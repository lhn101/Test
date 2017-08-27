function [estim]=online_calibration_NR_EM(data,model,option,varargin)
index.plot_estim=0;        % smother option plot (estim.smooth) in plot_estimations.m will not be 
                           % used in the online procedure | save_estimation_results.m 
filter_smoother=0;         % filter_smoother=1 Kalman filter and smoother are employed at the same time 
args = varargin;
for t=1:2:length(args)
    switch args{t}
        case 'filter_smoother', filter_smoother = args{t+1};
        otherwise, error(['unrecognized argument ' args{t}])
    end
end


%% Preallocating the memory of estimation results -> Reduce memory operations
[estim.filter]=preallocate_memory(data,model,option);
if filter_smoother==1
    [estim.smoother]=preallocate_memory(data,model,option);
end

%% Step 1 - Initialize the estimation and calibration with the first traning period
training_end_idx_ref= option.training_end_idx;% Identify the first end-point of training period 
                                              % in order to evaluate remaining time in Waitbar 
parameter_ref=model.parameter; 
[optim]=NR_EM_V2(data,model,option);
model.parameter=optim.parameter_opt;

data_train=data;
data_train.timestamps=data.timestamps(option.training_start_idx:option.training_end_idx);
data_train.values=data.values(option.training_start_idx:option.training_end_idx,:);
if isfield(data,'ref')% Reference data
    data_train.ref=data.ref(option.training_start_idx:option.training_end_idx,:);
end
data_train.dt_steps=data.dt_steps(option.training_start_idx:option.training_end_idx);
data_train.nb_steps=length(option.training_start_idx:option.training_end_idx);

estim_train.filter= state_estimation(data_train,model,option,'disp_flag',0);
if filter_smoother==1
    estim_train.smoother= state_estimation(data_train,model,option,'smooth',1,'disp_flag',0);
end

% Save estimation results of 1st training period
index.idx_1=option.training_start_idx;
index.idx_2=option.training_end_idx;
index.idx_3=index.idx_1;
index.idx_4=index.idx_2;
[estim.filter]=save_estimation_results(estim.filter,estim_train.filter,model,index);
if filter_smoother==1
    [estim.smoother]=save_estimation_results(estim.smoother,estim_train.smoother,model,index);
end
% Plot
option.prev_training_end_idx=option.training_end_idx;
plot_online_estimations(estim.filter,data,model,option,0)

%% Step 2 - Estimate hidden state and calibrate parameter recursively
loop=0;
h=waitbar(0,'Computing remaining time...','Name','States estimation progression');
time_tot=0;


while option.training_end_idx<length(data.timestamps)
    loop=loop+1;
    tic;
    
    %% Time and model setup
    [option]=time_selection(data,option);
    nb_data_input=option.training_end_idx-option.prev_training_end_idx-1;
    model=recursive_setup(model,estim.filter,option);
    
    %% Initialize the parameters involving to Switching Kalman Filter (SKF)
    idx_initParam=find(~strcmp(model.param_properties(:,3),'1')); % SKF's parameters to be initialized
    idx_initAR=find(contains(model.param_properties(:,2),'AR'));  % AR parameters to be initialized
    idx_initSigma=find(model.parameter<1E-15);                     
    idx_initParam=[idx_initParam;idx_initAR;idx_initSigma];
    model.parameter(idx_initParam)=parameter_ref(idx_initParam);  % Initialize the parameters related to abnormal model
    
    %% NR_EM optimization
    [optim]=NR_EM_V2(data,model,option);
    model.parameter=optim.parameter_opt;
    
    %% Resize the datased with respect to the choosen training period for state_estimation.m
    data_train=data;
    data_train.timestamps=data.timestamps(option.training_start_idx:option.training_end_idx);
    data_train.values=data.values(option.training_start_idx:option.training_end_idx,:);
    if isfield(data,'ref') % Reference data
        data_train.ref=data.ref(option.training_start_idx:option.training_end_idx,:);
    end
    data_train.dt_steps=data.dt_steps(option.training_start_idx:option.training_end_idx);
    data_train.nb_steps=length(option.training_start_idx:option.training_end_idx);
    estim_train.filter= state_estimation(data_train,model,option,'disp_flag',0);
    if filter_smoother==1
        estim_train.smoother= state_estimation(data_train,model,option,'smooth',1,'disp_flag',0);
    end
    
    %% Save estimation results
    index.idx_1=option.prev_training_end_idx+1;
    index.idx_2=option.training_end_idx;
    index.idx_3=data_train.nb_steps-nb_data_input;
    index.idx_4=data_train.nb_steps;
    [estim.filter]=save_estimation_results(estim.filter,estim_train.filter,model,index);
    if filter_smoother==1
        [estim.smoother]=save_estimation_results(estim.smoother,estim_train.smoother,model,index);
    end
    %% Waitbar update
    time_loop=toc;
    time_tot=time_tot+time_loop;
    time_rem=time_tot/(option.training_end_idx-training_end_idx_ref)*(length(data.timestamps)-option.training_end_idx)/60;
    
    %% Plot
    plot_online_estimations(estim.filter,data,model,option,loop)
    waitbar((option.training_end_idx-training_end_idx_ref)/(length(data.timestamps)-training_end_idx_ref),h,['Estimation loop #',sprintf('%0.0f',loop),':',...
            sprintf('%8.0f',(option.training_end_idx-training_end_idx_ref)),'/',sprintf('%2.0f',(length(data.timestamps)-training_end_idx_ref)),' ',...
            ' [',sprintf('%0.1f',time_rem),' ','min left]']);

end
close(h)
end

function [option]=time_selection(data,option)
datetime_vec=datevec(data.timestamps);
option.prev_training_end_idx=option.training_end_idx;     % Save the previous end-point training
option.prev_training_start_idx=option.training_start_idx; % Save the previous start-point training

option.training_end_idx=find(abs(data.timestamps-data.timestamps(option.prev_training_end_idx+1)-option.onl_estim_window+1)==...
                         min(abs(data.timestamps-data.timestamps(option.prev_training_end_idx+1)-option.onl_estim_window+1)),1,'first');
option.training_start_num=datenum(datetime(datetime_vec(option.training_end_idx,:)) - caldays(option.training_period(2)));
option.training_start_idx=find(abs(data.timestamps(1:option.training_end_idx)-option.training_start_num)==...
                           min(abs(data.timestamps(1:option.training_end_idx)-option.training_start_num)),1,'first');

if (option.training_end_idx-option.training_start_idx)< 0.75*option.nb_data_training % Make sure that there are enough data points for training period 
    option.training_start_idx= option.training_end_idx-0.75*option.nb_data_training;
    if option.training_start_idx<0
        error('Last data-point for training period is not properly extracted | time_selection.m')
    end
end    
end

function model=recursive_setup(model,estim,option)
% Hidden states & state probabilities
for i=1:model.nb_class
    model.initX{i}=estim.x_M{i}(:,option.training_start_idx-1); 
    model.initV{i}=estim.V_M{i}(:,:,option.training_start_idx-1);
    model.initS{i}=estim.S(option.training_start_idx-1,i);
end
% UD filter
if strcmp(option.method,'UD')
    model.U=cell(model.nb_class,model.nb_class,1);
    model.D=cell(model.nb_class,model.nb_class,1);
    for j=1:model.nb_class
        for i=1:model.nb_class
            model.U{i,j,1}=estim.U{i,j,option.training_start_idx};
            model.D{i,j,1}=estim.D{i,j,option.training_start_idx};
        end
    end
end
end

function [estim]=preallocate_memory(data,model,option)
ss=size(model.hidden_states_names{1},1);
estim.x=zeros(ss,data.nb_steps);
estim.V=zeros(ss,data.nb_steps);
estim.y=zeros(model.nb_obs,data.nb_steps);
estim.Vy=zeros(model.nb_obs,data.nb_steps);
estim.S=zeros(data.nb_steps,model.nb_class);
estim.parameter=zeros(data.nb_steps,size(model.parameter,1));
estim.log_lik=zeros(data.nb_steps,1);

if strcmp(option.method,'UD')
    estim.U=cell(model.nb_class,model.nb_class,data.nb_steps);
    estim.D=cell(model.nb_class,model.nb_class,data.nb_steps);
end
estim.x_M = cell(1,model.nb_class);
estim.V_M = cell(1,model.nb_class);
for j=1:model.nb_class
    estim.x_M{j} = zeros(ss,data.nb_steps);
    estim.V_M{j} = zeros(ss, ss,data.nb_steps);
end
end

function [estim]=save_estimation_results(estim,estim_train,model,index)

% Hidden states & state probilities
for j=1:model.nb_class
    estim.x_M{j}(:,index.idx_1:index.idx_2)=estim_train.x_M{j}(:,index.idx_3:index.idx_4);
    estim.V_M{j}(:,:,index.idx_1:index.idx_2)=estim_train.V_M{j}(:,:,index.idx_3:index.idx_4);
end
estim.x(:,index.idx_1:index.idx_2)=estim_train.x(:,index.idx_3:index.idx_4);
estim.V(:,index.idx_1:index.idx_2)=estim_train.V(:,index.idx_3:index.idx_4);
estim.y(:,index.idx_1:index.idx_2)=estim_train.y(:,index.idx_3:index.idx_4);
estim.Vy(:,index.idx_1:index.idx_2)=estim_train.Vy(:,index.idx_3:index.idx_4);
estim.S(index.idx_1:index.idx_2,:)=estim_train.S(index.idx_3:index.idx_4,:);

% UD filter
if isfield(estim,'U')
    estim.U(:,:,index.idx_1:index.idx_2)=estim_train.U(:,:,index.idx_3:index.idx_4);
    estim.D(:,:,index.idx_1:index.idx_2)=estim_train.D(:,:,index.idx_3:index.idx_4);
end
estim.smooth=index.plot_estim;

for t=index.idx_1:index.idx_2
    estim.parameter(t,:)=model.parameter;
    estim.log_lik(t,:)=estim_train.LL;
end
end