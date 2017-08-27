function [optim]=NR_EM_V2(data,model,option)
if ~isfield(option,'disp_flag')
    disp_flag=1;
end
%% Resize the dataset according to the choosen training period

data.timestamps=data.timestamps(option.training_start_idx:option.training_end_idx);
data.values=data.values(option.training_start_idx:option.training_end_idx,:);
data.dt_steps=data.dt_steps(option.training_start_idx:option.training_end_idx);
data.nb_steps=length(option.training_start_idx:option.training_end_idx);

%% Identify parameters to be optimized
parameter_search_idx=find(~all(isnan(reshape([model.param_properties{:,5}],2,size(model.param_properties,1))'),2));

%% Initialize transformed model parameters
parameter=model.parameter;
parameter_TR=model.parameter;
for i=1:length(parameter_search_idx)
    idx=parameter_search_idx(i);
    [fct_TR,~,~,~]=parameter_transformation_fct(model,idx);
    parameter_TR(idx)=fct_TR(parameter(idx));
end
parameter_TR_ref=parameter_TR;

%% Analysis parameters
nb_levels_lambda_ref=4;
delta_grad= 1E-3*ones(size(model.parameter));
convergence_tolerance=1E-7;

if disp_flag==1
    %% Diaplay analysis parameters
    disp('----------------------------------------------------------------------------------------------')
    disp('Estimate the parameters for the Bayesian Dynamic Linear Model')
    disp('----------------------------------------------------------------------------------------------')
    disp(' ')
    disp('    \\start Newton-Raphson maximization algorithm (finite difference method)')
    disp(' ')
    disp(['      Training period:                                       ' num2str(option.training_period(1)) '-' num2str(option.training_period(2)) ' [days]'])
    disp(['      Maximal number of iteration:                           ' num2str(option.iteration_limit_calibration)])
    disp(['      Total time limit for calibration:                      ' num2str(option.time_limit_calibration) ' [min]'])
    disp(['      Convergence criterion:                                 ' num2str(convergence_tolerance) '*LL'])
    disp(['      Nb. of search levels for \lambda:                      ' num2str(nb_levels_lambda_ref) '*2'])
    disp(['      \delta for numerical gradient:                         ' num2str(delta_grad(1)) '*param'])
    disp(' ')
    disp('    ...in progress')
    disp(' ')
end

%% Matrices & parameter initilization
hessian_calculation_fail=zeros(1,numel(parameter_search_idx));
optimization_fail=zeros(1,numel(parameter_search_idx));
converged=zeros(1,numel(parameter_search_idx));

dll=1E6*ones(1,numel(parameter_search_idx));

grad_p_TR=zeros(size(parameter));
grad_p=zeros(size(model.parameter));
hessian_p_TR=zeros(size(parameter));
std_p=zeros(size(parameter));
std_p_TR=parameter_TR;

%% Log-likelihood initialization
[~,~,~,~,log_lik_0,~,~]=SKF(data,model,option);
if isinf(log_lik_0)
    disp('warning LL0=-inf | NR_EM.m')
end
if disp_flag==1
    disp(['           Initial LL: ' num2str(log_lik_0)])
end
name_idx_1='';
name_idx_2='';
for i=parameter_search_idx'
    name_p1{i}=[model.param_properties{i,1}];
    if ~isempty(model.param_properties{i,4})
        temp=model.param_properties{i,4}(1);
    else
        temp='';
    end
    name_p2{i}=[model.param_properties{i,2}, '|M', model.param_properties{i,3},'|',temp];
    name_idx_1=[name_idx_1  name_p1{i} repmat(' ',[1,15-length(name_p1{i})]) ' '];
    name_idx_2=[name_idx_2  name_p2{i} repmat(' ',[1,15-length(name_p2{i})]) ' '];  
end
if disp_flag==1
    disp(['                       ' name_idx_2])
    disp(['      parameter names: ' name_idx_1])
    disp(sprintf(['       initial values: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
end
param_idx_prev=1;

%% NR Optimization loops
tic; % time counter initialization
time_loop=0;
search_loop=1;
while and(search_loop<=option.iteration_limit_calibration, time_loop<(option.time_limit_calibration*60))
    if disp_flag==1
        disp(' ')
    end
    parameter_ref=parameter;
    
    %% Select the parameter which previously led to the highest change in LL
    rand_sample=rand;
    dll_cumsum=cumsum(dll)/sum(dll);
    dll_rank=dll_cumsum-rand_sample;
    dll_rank(dll_rank<0)=inf;
    dll_rank(hessian_calculation_fail>10)=Inf; % Remove the parameters that have small impact on log-likelihood computation
    if all(dll_rank==Inf)
        param_idx=find(~converged,1,'first');
    else
        param_idx=find((dll_cumsum-rand_sample)==min(dll_rank),1,'first');
    end   
    param_idx_loop=parameter_search_idx(param_idx);
%     idx_initParam=parameter_search_idx(abs(parameter_ref(parameter_search_idx))<1E-12);
%     parameter(idx_initParam)=1E-8;           % <1E-15 the computational accurancy will be an issue
%                                               % Also allow speeding up the EM_NR algorithm
%     parameter_ref(idx_initParam)=1E-8;
    param=parameter_ref(param_idx_loop);
    param_TR=parameter_TR(param_idx_loop);
    [fct_TR,fct_inv_TR,grad_OR2TR,hessian_OR2TR]=parameter_transformation_fct(model,param_idx_loop);
    if disp_flag==1
        disp('--------------------------')
        disp(['    Loop #' num2str(search_loop) ' : ' [name_p2{param_idx_loop} '|' name_p1{param_idx_loop}]])
    end
    H_test=1;
    delta_grad_loop=delta_grad(param_idx_loop)*abs(param);
    loop_count=1;
    while H_test                              % loop until the hessian can be calculated
        loop_count=loop_count +1;
        skip_loop_H=0;
        if loop_count >15
            if disp_flag==1
                disp('           Warning:>15 failed attempts to compute the Hessian')
            end
            skip_loop_H=1;
            hessian_calculation_fail(param_idx)=hessian_calculation_fail(param_idx)+1;
            parameter(param_idx_loop)=parameter_ref(param_idx_loop);
            parameter_TR(param_idx_loop)=parameter_TR_ref(param_idx_loop);
            if hessian_calculation_fail(param_idx)>3
                if disp_flag==1
                    disp('           Warning: 3nd discontinued fails to compute the Hessian  -> converged=1')
                end
                converged(param_idx)=1;
            end
            dll(param_idx)=abs(convergence_tolerance*log_lik_0);
            break           
        end
        
        %% Boundary condition validity in original space
        min_p=model.param_properties{param_idx_loop,5}(1);
        max_p=model.param_properties{param_idx_loop,5}(2);
        delta_grad_OR=delta_grad_loop;
        delta_test=1;
        delta_loop=0;
        while delta_test
            delta_loop=delta_loop+1;
            if delta_loop>10
                disp('           Warning: delta test has failed')
                skip_loop_H=1;
                break
            end
            if or((param-delta_grad_OR)<min_p,(param+delta_grad_OR>max_p))
                delta_grad_OR=delta_grad_OR/2;
                delta_grad(param_idx_loop)=delta_grad_OR/abs(param);
                delta_test=1;
                if disp_flag==1
                    disp('           Warning: Parameter across its boundary -> delta_grad_OR/2')
                end
            else
                delta_test=0;
            end
        end
        grad_TR=@(p) (fct_TR(p)*delta_grad_loop)/(fct_inv_TR((fct_TR(p)*(1+delta_grad_loop)))-p);
        grad_inv_TR=@(p) 1/grad_TR(p);
        
        if skip_loop_H==1
            hessian_calculation_fail(param_idx)=hessian_calculation_fail(param_idx)+1;
            parameter(param_idx_loop)=parameter_ref(param_idx_loop);
            parameter_TR(param_idx_loop)=parameter_TR_ref(param_idx_loop);
            if hessian_calculation_fail(param_idx)>3
                if disp_flag==1
                    disp('           Warning: delta test has failed 3 times  -> converged=1')
                end
                converged(param_idx)=1;
            end
            dll(param_idx)=abs(convergence_tolerance*log_lik_0);
            break
        end
        %% Log-likehood calculation
        log_lik_s=[0,0];
        % Parallele computing
        if strcmp(option.computing_method,'parallel')            
            p_LL=[model.parameter,model.parameter];
            p_LL(param_idx_loop,1)=param-delta_grad_OR;
            p_LL(param_idx_loop,2)=param+delta_grad_OR;
            s=cell(1,2);
            m_1=model;
            m_2=model;
            m_1.parameter=p_LL(:,1);
            m_2.parameter=p_LL(:,2);
            s{1}=m_1;
            s{2}=m_2;
            parfor i=1:2
                [~,~,~,~,log_lik_s(i),~,~]=SKF(data,s{i},option); % LL calcultation
            end
        % Serial computing
        elseif strcmp(option.computing_method,'serial')
            model.parameter(param_idx_loop)=param-delta_grad_OR;
            [~,~,~,~,log_lik_s(1),~,~]=SKF(data,model,option);
            model.parameter(param_idx_loop)=param+delta_grad_OR;
            [~,~,~,~,log_lik_s(2),~,~]=SKF(data,model,option);           
        end
        log_lik_1=log_lik_s(1);
        log_lik_2=log_lik_s(2);
        
        if any(~isreal([log_lik_0,log_lik_1,log_lik_2]))
            disp('           Warning: LL is complex -> LL=real(LL)')
            log_lik_0=real(log_lik_0);
            log_lik_1=real(log_lik_1);
            log_lik_2=real(log_lik_2);
        end
        %% Size-step optimization ensure that the gradient can not change too rapidly or slowly
        if abs((log_lik_2-log_lik_1))<=(1E-2*delta_grad_OR*2)  %Need to identify the compatible threshold  
            delta_grad_loop=2*delta_grad_loop;
            delta_grad(param_idx_loop)=2*delta_grad(param_idx_loop);
            if disp_flag==1
                disp(['           Warning: Step size has increased -> delta_grad_loop: ' sprintf('%#-+8.2e',delta_grad_loop) ' (delta_grad_loop=2*delta_grad_loop)'])
            end
        elseif abs((log_lik_2-log_lik_1)/log_lik_0)>=1E-2 % Need to identify the compatible threshold
            delta_grad_loop=delta_grad_loop/2;
            delta_grad(param_idx_loop)=delta_grad(param_idx_loop)/2;
            if disp_flag==1
                disp(['           Warning: Step size has decreased -> delta_grad_loop: ' sprintf('%#-+8.2e',delta_grad_loop) ' (delta_grad_loop=delta_grad_loop/2)'])
            end
        else
        %% Grad & hessian calculation
            grad_p_TR(param_idx_loop)=((log_lik_2-log_lik_0)/delta_grad_OR)*grad_OR2TR(param_TR);
            hessian_p_TR_loop=((log_lik_2-2*log_lik_0+log_lik_1)/(delta_grad_OR^2))*grad_OR2TR(param_TR)^2 ...
            +((log_lik_2-log_lik_0)/delta_grad_OR)*hessian_OR2TR(param_TR);
        
            if hessian_calculation_fail(param_idx)>0
                hessian_calculation_fail(param_idx)=0;
            end 
            if grad_p_TR(param_idx_loop)==0
                if disp_flag==1
                    disp('           Warning: grad(p)=0 -> delta_grad_loop=delta_grad_loop*10^n') %gradient cannot be calculated
                end
                delta_grad_loop=delta_grad_loop*10^loop_count;  
            elseif hessian_p_TR_loop<0 
                hessian_p_TR(param_idx_loop)=hessian_p_TR_loop;
                delta_p_TR=-grad_p_TR(param_idx_loop)/hessian_p_TR(param_idx_loop);% Standard Newton-Raphson
                H_test=0;
            elseif hessian_p_TR_loop>0 
                hessian_p_TR(param_idx_loop)=hessian_p_TR_loop;
                if model.param_properties{param_idx_loop,5}(1)==0&&model.param_properties{param_idx_loop,5}(2)==inf&&param<1E-7
                    delta_p_TR=sign(grad_p_TR(param_idx_loop))*0.25*param_TR;
                else
                    delta_p_TR=grad_p_TR(param_idx_loop)/hessian_p_TR(param_idx_loop);%Gradient descent must be used to reach the maxima
                end
                if disp_flag==1
                    disp('           Warning: Hessian is positive ')
                end
                H_test=0; 
            elseif hessian_p_TR_loop==0 
                if disp_flag==1
                    disp('           Warning: H(p)=0 -> delta_grad_loop=delta_grad_loop*2') %Hessian cannot be calculated
                end
                delta_grad_loop=delta_grad_loop*2;                
            end    
        end
    end
    if skip_loop_H~=1
        std_p_TR_loop=sqrt(-1/hessian_p_TR(param_idx_loop));
        std_p_loop=abs(grad_inv_TR(param))*std_p_TR_loop;%Linearized Laplace approximation of parameter standard deviation
        grad_p(param_idx_loop)=abs(grad_inv_TR(param))*grad_p_TR(param_idx_loop);
        if isreal(std_p_loop)
            if std_p_loop>0
                std_p(param_idx_loop)=std_p_loop;
                std_p_TR(param_idx_loop)=std_p_TR_loop;
            else
                std_p(param_idx_loop)=nan;
            end
        else
            std_p(param_idx_loop)=nan;
        end
        
        delta_p=fct_inv_TR(param_TR+delta_p_TR)-param;
        if disp_flag==1
            disp(['       delta_param: ' num2str(delta_p)])
        end
        try
        %% LL test
            parameter_TR(param_idx_loop)=param_TR+delta_p_TR;
            parameter(param_idx_loop)=fct_inv_TR(parameter_TR(param_idx_loop));
            model.parameter=parameter;
            [~,~,~,~,log_lik_test,~,~]=SKF(data,model,option);
        catch err
            log_lik_test=-inf;
        end
        loop_converged=1;
        delta_p_TR_ref=delta_p_TR;
    else
        loop_converged=0; % Optimization has failed the actual parameter -> move to next parameter
    end
    %% Check the validity of delta_p
    n=0;
    reverse=0;
    n_ref=0;
    nb_levels_lambda=nb_levels_lambda_ref;
    while loop_converged 
        n=n+1;
        if log_lik_test>log_lik_0
            loop_converged=0;
            converged(param_idx)=and(~isnan(std_p(param_idx_loop)),abs(log_lik_test-log_lik_0)<abs(convergence_tolerance*log_lik_0));
            dll(param_idx)=log_lik_test-log_lik_0; % Record the change in the log-likelihood
            if dll(param_idx)<abs(convergence_tolerance*log_lik_0)
                if isnan(std_p(param_idx_loop))
                    dll(param_idx)=10*abs(convergence_tolerance*log_lik_0);
                else
                    dll(param_idx)=abs(convergence_tolerance*log_lik_0);
                end
            end
            log_lik_0=log_lik_test;
            parameter_TR_ref(param_idx_loop)=parameter_TR(param_idx_loop);
            parameter_ref(param_idx_loop)=parameter(param_idx_loop);
            if disp_flag==1
                disp(['    log-likelihood: ' num2str(log_lik_0)])
                disp(['      param change: ' num2str(param) ' -> ' num2str(parameter(param_idx_loop))])
                disp(' ')
                disp(['                    ' name_idx_2])
                disp(['   parameter names: ' name_idx_1])
                disp(sprintf(['    current values: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
                disp(sprintf(['  current f.o. std: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],std_p(parameter_search_idx)))
                disp(sprintf(['      previous dLL: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],dll))
                disp(sprintf(['          gradient: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],grad_p(parameter_search_idx)))
                disp(sprintf(['           hessian: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],hessian_p_TR(parameter_search_idx)))
                disp(sprintf(['         converged: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],converged))
                disp(sprintf(['      hessian_fail: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],hessian_calculation_fail))
            end
            optimization_fail=optimization_fail*0;
            
            if converged(param_idx)
                break
            else
                converged(param_idx)=0;               
            end
            
        else
            converged(param_idx)=0;
            if n>nb_levels_lambda
                if disp_flag==1
                    disp(' ')
                    disp('      ...optimization loop has failed')
                end
                parameter(param_idx_loop)=parameter_ref(param_idx_loop);
                parameter_TR(param_idx_loop)=parameter_TR_ref(param_idx_loop);
                optimization_fail(param_idx)=optimization_fail(param_idx)+1;
                if optimization_fail(param_idx)>3
                    converged(param_idx)=1;
                    dll(param_idx)=abs(convergence_tolerance*log_lik_0);
                else
                    dll(param_idx)=2*abs(convergence_tolerance*log_lik_0);
                end
                break
            elseif n<nb_levels_lambda
                delta_p_TR=delta_p_TR/(2^n);
                delta_p=fct_inv_TR(param_TR+delta_p_TR)-param;
                if disp_flag==1
                    disp(['           Warning: log-likelihood has decreased -> delta_param: ' sprintf('%#-+8.2e',delta_p) ' (delta_p_TR=delta_p_TR/(2^' num2str(n) '))'])
                end
            elseif n==nb_levels_lambda
                if reverse==0
                    reverse=1;
                    n=n_ref;
                    delta_p_TR=-delta_p_TR_ref;
                    delta_p=fct_inv_TR(param_TR+delta_p_TR)-param;
                elseif reverse==1
                    reverse=0;
                    n_ref=n;
                    nb_levels_lambda=2*nb_levels_lambda_ref;
                    delta_p_TR=delta_p_TR_ref;
                    delta_p=fct_inv_TR(param_TR+delta_p_TR)-param;
                end
                if disp_flag==1
                    disp(['           Warning: log-likelihood has decreased -> delta_param: ' sprintf('%#-+8.2e',delta_p) ' (delta_p_TR=-delta_p_TR)'])
                end
                
            end
            try
                %% LL test
                parameter_TR(param_idx_loop)=param_TR+delta_p_TR;
                parameter(param_idx_loop)=fct_inv_TR(parameter_TR(param_idx_loop));
                model.parameter=parameter;
                [~,~,~,~,log_lik_test,~,~]=SKF(data,model,option);
            catch err
                log_lik_test=-inf;
            end
        end
    end
    search_loop=search_loop+1;
    if all(optimization_fail>0)
        if disp_flag==1
            disp(' ')
            disp('          WARNING: the optimization has failed for every parameter')
        end
        break
    elseif all(converged)
        if disp_flag==1
            disp(' ')
            disp('             DONE: the optimization has converged for all parameters')
        end
        break
    end
    time_loop=toc;   
end
if disp_flag==1
    if ~or(all(optimization_fail),all(converged))
        disp(' ')
        disp(' ')
        disp(' ')
        if time_loop>(option.time_limit_calibration*60)
            disp(['          WARNING : the optimization has reached the maximum allowed time (' num2str(option.time_limit_calibration) ' [min]) without convergence'])
        else
            disp(['          WARNING : the optimization has reached the maximum number of loops (' num2str(option.iteration_limit_calibration) ') without convergence'])
        end
    end
end

if disp_flag==1
    %% Display final results
    disp(' ')
    disp(' ----------------------')
    disp('    Final results')
    disp(' ----------------------')
    disp(['     log-likelihood: ' num2str(log_lik_0)])
    disp(['                     ' name_idx_2])
    disp(['    parameter names: ' name_idx_1])
    disp(sprintf(['   optimized values: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
    disp(sprintf(['   f.o. Laplace std: ' repmat(['%#-+15.2e' ' '],[1,length(parameter)])],std_p(parameter_search_idx)))
    disp(' ')
    disp('***************************************** End *****************************************')
    disp(' ')
    disp(' ')
end
optim.parameter_opt=parameter;
optim.std_parameter_opt=std_p;
optim.converged=converged;
optim.search_loop=search_loop;
optim.log_lik=log_lik_0;
end

function [fct_TR,fct_inv_TR,grad_OR2TR,hessian_OR2TR]=parameter_transformation_fct(model,param_idx_loop)
sigmoid=@(p) 1./(1+exp(-p));
sigmoid_inv=@(p) -log((1./p)-1);

p_TR_0=@(p) p;
p_TR_1=@(p) log(p);
p_TR_2=@(p,min_p,range_p) sigmoid_inv((p-min_p)/range_p);


p_inv_TR_0=@(p) p;
p_inv_TR_1=@(p) exp(p);
p_inv_TR_2=@(p,min_p,range_p) sigmoid(p)*range_p+min_p;


if model.param_properties{param_idx_loop,5}(1)==-inf&&model.param_properties{param_idx_loop,5}(2)==inf
    fct_TR=@(p) p_TR_0(p);
    fct_inv_TR=@(p) p_inv_TR_0(p);
    grad_OR2TR=@(p) 1;
    hessian_OR2TR=@(p) 1;
elseif model.param_properties{param_idx_loop,5}(1)==0&&model.param_properties{param_idx_loop,5}(2)==inf
    fct_TR=@(p) p_TR_1(p);
    fct_inv_TR=@(p) p_inv_TR_1(p);
    grad_OR2TR=@(p) exp(p);
    hessian_OR2TR=@(p) exp(p);
elseif isfinite(model.param_properties{param_idx_loop,5}(1))&&isfinite(model.param_properties{param_idx_loop,5}(2))
    range_p=model.param_properties{param_idx_loop,5}(2)-model.param_properties{param_idx_loop,5}(1);
    min_p=model.param_properties{param_idx_loop,5}(1);
    fct_TR=@(p) p_TR_2(p,min_p,range_p);
    fct_inv_TR=@(p) p_inv_TR_2(p,min_p,range_p);
    grad_OR2TR=@(p) sigmoid(p)*(1-sigmoid(p))*range_p;
    hessian_OR2TR=@(p) sigmoid(p)*(1-sigmoid(p))*(1-2*sigmoid(p))*range_p;
else
    error('Parameter bounds are not properly defined in: model.param_properties{:,5}')
end
end