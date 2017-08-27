function [optim]=NR_EM(data,model,option)
if ~isfield(option,'disp_flag')
    disp_flag=1;
end

%% Resize the datased according to the choosen training period.
data.timestamps=data.timestamps(option.training_start_idx:option.training_end_idx);
data.values=data.values(option.training_start_idx:option.training_end_idx,:);
data.dt_steps=data.dt_steps(option.training_start_idx:option.training_end_idx);
data.nb_steps=length(option.training_start_idx:option.training_end_idx);

%% identify parameters to be optimized
parameter_search_idx=find(~all(isnan(reshape([model.param_properties{:,5}],2,size(model.param_properties,1))'),2));

%% initialize transformed model parameters
parameter=model.parameter;
parameter_TR=model.parameter;
for i=1:length(parameter_search_idx)
    idx=parameter_search_idx(i);
    [fct_TR,~]=parameter_transformation_fct(model,idx);
    parameter_TR(idx)=fct_TR(parameter(idx));
end
parameter_TR_ref=parameter_TR;


%% Analysis parameters
nb_levels_lambda_ref=3;
delta_grad=1E-3*ones(1,length(parameter_search_idx));
convengerce_tolerance=1E-7;

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
    disp(['      Convergence criterion:                                 ' num2str(convengerce_tolerance) '*LL'])
    disp(['      Nb. of search levels for \lambda:                      ' num2str(nb_levels_lambda_ref) '*2'])
    disp(['      \delta for numerical gradient:                         ' num2str(delta_grad(1)) '*param'])
    disp(' ')
    disp('    ...in progress')
    disp(' ')
end

%% Matrices & parameter initialization
optimization_fail=zeros(1,numel(parameter_search_idx));
converged=zeros(1,numel(parameter_search_idx));

dll=1E6*ones(1,numel(parameter_search_idx));

grad_p_TR=zeros(size(parameter));
hessian_p_TR=zeros(size(parameter));
std_p=zeros(size(parameter));
std_p_TR=parameter_TR;
sign_delta_grad_loop=ones(size(parameter));

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
    name_p1{i}=[model.param_properties{i,1} ];
    if ~isempty(model.param_properties{i,4})
        temp=model.param_properties{i,4}(1);
    else
        temp=' ';
    end
    name_p2{i}=[model.param_properties{i,2},'|M',model.param_properties{i,3},'|',temp];
    
    name_idx_1=[name_idx_1  name_p1{i} repmat(' ',[1,13-length(name_p1{i})]) ' '];
    name_idx_2=[name_idx_2  name_p2{i} repmat(' ',[1,13-length(name_p2{i})]) ' '];
    
end
if disp_flag==1
    disp(['                       ' name_idx_2])
    disp(['      parameter names: ' name_idx_1])
    disp(sprintf(['       initial values: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
end
param_idx_prev=1;

%% NR Optimization loops
tic; %time counter initialization
time_loop=0;
search_loop=0;
while and(search_loop<=option.iteration_limit_calibration,time_loop<(option.time_limit_calibration*60))
    if disp_flag==1
        disp(' ')
    end
    parameter_ref=parameter;
    
    %% Select the parameter who previously led to the highest change in LL
    rand_sample=rand;
    dll_cumsum=cumsum(dll)/sum(dll);
    dll_rank=dll_cumsum-rand_sample;
    dll_rank(dll_rank<0)=inf;
    param_idx=find((dll_cumsum-rand_sample)==min(dll_rank),1,'first');
    
    param_idx_loop=parameter_search_idx(param_idx);
    param=parameter_ref(param_idx_loop);
    param_TR=parameter_TR(param_idx_loop);
    [fct_TR,fct_inv_TR]=parameter_transformation_fct(model,param_idx_loop);
    if disp_flag==1
        disp(' ----------------------')
        disp(['    Loop #' num2str(search_loop) ' : ' [name_p2{param_idx_loop} '|' name_p1{param_idx_loop}]])
        disp(' ')
    end
    
    H_test=1;
    delta_grad_loop=sign_delta_grad_loop(param_idx_loop)*delta_grad(param_idx)*abs(param_TR);
    if delta_grad_loop==0
        delta_grad_loop=delta_grad(param_idx);
    end
    loop_count=0;
    while H_test %loop until the hessian can be calculated
        loop_count=loop_count+1;
        if loop_count>5;
            disp(' >5 failed attempts to compute the Hessian -> \delta_p=p')
            hessian_p_TR(param_idx_loop)=1;
            delta_p_TR=param_TR;
            break
        end
        
        delta_grad_TR=delta_grad_loop;
        grad_TR=@(p) (fct_TR(p)*delta_grad_loop)/(fct_inv_TR((fct_TR(p)*(1+delta_grad_loop)))-p);
        grad_inv_TR=@(p) 1/grad_TR(p);
        
        %% Parrallel Finite diff grad calculation
        log_lik_s=zeros(1,2);
        save_model=cell(1,2);
        p_LL=[model.parameter,model.parameter];
        p_LL(param_idx_loop,1)=fct_inv_TR(param_TR-delta_grad_TR);
        p_LL(param_idx_loop,2)=fct_inv_TR(param_TR+delta_grad_TR);
        m_1=model;
        m_2=model;
        m_1.parameter=p_LL(:,1);
        m_2.parameter=p_LL(:,2);
        save_model{1}=m_1;
        save_model{2}=m_2;
        parfor i=1:2
            [~,~,~,~,log_lik_s(i),~,~]=SKF(data,save_model{i},option); % LL calcultation
        end       
        log_lik_1=log_lik_s(1);
        log_lik_2=log_lik_s(2);
        
        if any(~isreal([log_lik_0,log_lik_1,log_lik_2]))
            disp('           Warning: LL is complex -> LL=real(LL)')
            log_lik_0=real(log_lik_0);
            log_lik_1=real(log_lik_1);
            log_lik_2=real(log_lik_2);
        end
        
        
        %% Grad & hessian calculation
        grad_p_TR(param_idx_loop)=(log_lik_1-log_lik_0)/delta_grad_TR;
        hessian_p_TR_loop=(log_lik_2-2*log_lik_0+log_lik_1)/(delta_grad_TR^2);
        
        if grad_p_TR(param_idx_loop)~=0
            if abs(grad_p_TR(param_idx_loop)*delta_grad(param_idx))<abs(log_lik_0*convengerce_tolerance)
                delta_grad(param_idx)=10*delta_grad(param_idx);
                disp('                    \delta grad(p)=10*\delta grad(p)') %gradient cannot be calculated
                disp(['                    \delta grad=[' num2str(delta_grad) ']*param'])
            elseif abs(grad_p_TR(param_idx_loop)*delta_grad(param_idx))>abs(convengerce_tolerance*1E4*log_lik_0)
                delta_grad(param_idx)=delta_grad(param_idx)/10;
                disp('                    \delta grad(p)=\delta grad(p)/10') %gradient cannot be calculated
            end
        end
        
        if grad_p_TR(param_idx_loop)==0
            if disp_flag==1
                disp('                   grad(p)=0 -> delta_grad_loop=delta_grad_loop*10^n') %gradient cannot be calculated
            end
            delta_grad_loop=delta_grad_loop*10^loop_count;
        elseif hessian_p_TR_loop<0
            hessian_p_TR(param_idx_loop)=hessian_p_TR_loop;
            delta_p_TR=-grad_p_TR(param_idx_loop)/hessian_p_TR(param_idx_loop);  %Standard Newton-Raphson
            H_test=0;
        elseif hessian_p_TR_loop>0
            delta_p_TR=grad_p_TR(param_idx_loop)/hessian_p_TR(param_idx_loop);  %Newton-Raphson in reverse direction
            H_test=0;
            if disp_flag==1
                disp('                    saddle point H(p)>0 -> delta_p_TR=-delta_p_TR')
            end
        elseif hessian_p_TR_loop==0
            if disp_flag==1
                disp('                    H(p)=0 -> delta_grad_loop=delta_grad_loop*2') %Hessian cannot be calculated
            end
            delta_grad_loop=delta_grad_loop*2;
            delta_p_TR=0;
        else
            delta_p_TR=0;
        end
    end
    std_p_TR_loop=sqrt(-1/hessian_p_TR(param_idx_loop));
    std_p_loop=abs(grad_inv_TR(param))*std_p_TR_loop; %Linearized Laplace approximation of parameter standard deviation
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
    
    if abs(param_TR+delta_p_TR)>20 % avoid steps that leads to a drastic change in parameter value
        delta_p_TR=sign(delta_p_TR)*0.5*(20-abs(param_TR));
        if disp_flag==1
            disp(['                     abs(param_TR+delta_p_TR)>20 -> delta_p_TR=sign(delta_p_TR)*20 '])
        end
    end
    
    delta_p=param-fct_inv_TR(param_TR+delta_p_TR);
    if delta_p==0
        delta_p=abs(param_TR)*sign_delta_grad_loop(param_idx_loop)*delta_grad;
        if disp_flag==1
            disp('                     delta_p=0 -> delta_p=grad_inv_TR*grad_p_TR')
        end
    end
    
    if abs(delta_p)>1*abs(param) % avoid steps that leads to a drastic change in parameter value
        delta_p=sign(delta_p)*1*abs(param);
        if disp_flag==1
            disp(['                     delta_param > abs(param) -> delta_p=sign(delta_p)*abs(param) '])
        end
    end
    
    if abs(param+delta_p)<0.01*abs(param) % avoid steps that leads to a zero parameter value
        delta_p=0.99*delta_p;
        if disp_flag==1
            disp(['                     abs(param+delta_p)<0.01*abs(param) -> delta_p=0.99*delta_p '])
        end
    end
    if disp_flag==1
        disp(['       delta_param: ' num2str(delta_p)])
    end
    try
        %% LL test
        parameter(param_idx_loop)=param+delta_p;
        parameter_TR(param_idx_loop)=fct_TR(parameter(param_idx_loop));
        if isinf(parameter_TR(param_idx_loop))
            parameter_TR(param_idx_loop)=sign(parameter_TR(param_idx_loop))*1E6;
            sign_delta_grad_loop(param_idx_loop)=-sign(parameter_TR(param_idx_loop));
        end
        model.parameter=parameter;
        [~,~,~,~,log_lik_test,~,~]=SKF(data,model,option);
    catch err
        log_lik_test=-inf;
    end
    
    %% Check the validity of delta_p
    loop_converged=1;
    n=0;
    reverse=0;
    n_ref=0;
    delta_p_TR_ref=delta_p_TR;
    nb_levels_lambda=nb_levels_lambda_ref;
    while loop_converged
        n=n+1;
        if log_lik_test>log_lik_0
            loop_converged=0;
            converged(param_idx)=and(~isnan(std_p(param_idx_loop)),abs(log_lik_test-log_lik_0)<abs(convengerce_tolerance*log_lik_0));
            dll(param_idx)=log_lik_test-log_lik_0; % Record the change in the log-likelihood
            if dll(param_idx)<abs(convengerce_tolerance*log_lik_0)
                if isnan(std_p(param_idx_loop))
                    dll(param_idx)=10*abs(convengerce_tolerance*log_lik_0);
                else
                    dll(param_idx)=abs(convengerce_tolerance*log_lik_0);
                end
            end
            log_lik_0=log_lik_test;
            parameter_TR_ref(param_idx_loop)=parameter_TR(param_idx_loop);
            if disp_flag==1
                disp(['    log-likelihood: ' num2str(log_lik_0)])
                disp(['      param change: ' num2str(param) ' -> ' num2str(parameter(param_idx_loop))])
                disp(' ')
                disp(['                    ' name_idx_2])
                disp(['   parameter names: ' name_idx_1])
                disp(sprintf(['    current values: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
                disp(sprintf(['  current f.o. std: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],std_p(parameter_search_idx)))
                disp(sprintf(['      previous dLL: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],dll))
            end
            optimization_fail(param_idx)=optimization_fail(param_idx)*0;
            if converged(param_idx)
                break
            end
        else
            converged(param_idx)=0;
            if n>nb_levels_lambda
                if disp_flag==1
                    disp(' ')
                    disp(['    ...optimization loop has failed'])
                end
                parameter(param_idx_loop)=parameter_ref(param_idx_loop);
                parameter_TR(param_idx_loop)=parameter_TR_ref(param_idx_loop);
                optimization_fail(param_idx)=optimization_fail(param_idx)+1;
                if optimization_fail(param_idx)>5
                    converged(param_idx)=1;
                    dll(param_idx)=abs(convengerce_tolerance*log_lik_0);
                else
                    dll(param_idx)=2*abs(convengerce_tolerance*log_lik_0);
                end
                break
                
            elseif n<nb_levels_lambda
                if n==1
                    delta_p_ref=delta_p;
                end
                delta_p=delta_p_ref/10^n;
                if disp_flag==1
                    disp(['                     log-likelihood has decreased -> delta_param: ' sprintf('%#-+8.2e',delta_p) ' (delta_p=delta_p/(10^' num2str(n) '))'])
                end
                
            elseif n==nb_levels_lambda&&reverse==0
                reverse=1;
                n=0;
                delta_p=-delta_p_ref;
                if disp_flag==1
                    disp(['                     log-likelihood has decreased -> delta_param: ' sprintf('%#-+8.2e',delta_p) ' (delta_p=-delta_p)'])
                end
                
            end
            try
                %% LL test
                parameter(param_idx_loop)=param+delta_p;
                parameter_TR(param_idx_loop)=param_TR+delta_p_TR;
                if isinf(parameter_TR(param_idx_loop))
                    parameter_TR(param_idx_loop)=sign(parameter_TR(param_idx_loop))*1E6;
                    sign_delta_grad_loop(param_idx_loop)=-sign(parameter_TR(param_idx_loop));
                end
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
%% Save parameters, hessian and std(p)
std_p=std_p;
hessian_p_TR=hessian_p_TR;
if ~isreal(sqrt(-1./hessian_p_TR(parameter_search_idx)))
    disp('img hessian')
end
dll=dll;
parameter_search_idx=parameter_search_idx;
name_idx_1=name_idx_1;

if disp_flag==1
    %% Display final results
    disp(' ')
    disp(' ----------------------')
    disp('    Final results')
    disp(' ----------------------')
    disp(['     log-likelihood: ' num2str(log_lik_0)])
    disp(['                     ' name_idx_2])
    disp(['    parameter names: ' name_idx_1])
    disp(sprintf(['   optimized values: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],parameter(parameter_search_idx)))
    disp(sprintf(['   f.o. Laplace std: ' repmat(['%#-+13.2e' ' '],[1,length(parameter)])],std_p(parameter_search_idx)))
    disp(' ')
end
optim.parameter_opt=parameter;
optim.std_parameter_opt=std_p;
optim.converged=converged;
optim.search_loop=search_loop;
optim.log_lik=log_lik_0;
end

function [fct_TR,fct_inv_TR]=parameter_transformation_fct(model,param_idx_loop)
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
elseif model.param_properties{param_idx_loop,5}(1)==0&&model.param_properties{param_idx_loop,5}(2)==inf
    fct_TR=@(p) p_TR_1(p);
    fct_inv_TR=@(p) p_inv_TR_1(p);
elseif isfinite(model.param_properties{param_idx_loop,5}(1))&&isfinite(model.param_properties{param_idx_loop,5}(2))
    range_p=model.param_properties{param_idx_loop,5}(2)-model.param_properties{param_idx_loop,5}(1);
    min_p=model.param_properties{param_idx_loop,5}(1);
    fct_TR=@(p) p_TR_2(p,min_p,range_p);
    fct_inv_TR=@(p) p_inv_TR_2(p,min_p,range_p);
else
    error('Parameter bounds are not properly defined in: model.param_properties{:,5}')
end
end