function [model,data,option]=model_initialisation(model,data,option)
model.nb_class=size(model.components.block,2);        %Number of model classes
model.nb_obs=size(model.components.block{1},2);       %Number of observations

data.dt_mean=mean(diff(data.timestamps));        %Mean time step length
data.dt_steps=diff(data.timestamps);
unique_dt_steps=unique(data.dt_steps);
counts_dt_steps= [unique_dt_steps,histc(data.dt_steps(:),unique_dt_steps)];
data.dt_ref=counts_dt_steps(find(counts_dt_steps(:,2)==max(counts_dt_steps(:,2)),1,'first'),1);  %Define the reference time step as the most frequent
data.dt_steps=[data.dt_ref;data.dt_steps];       %Define time step vector
data.nb_steps=length(data.timestamps);           %Store the number of time steps

%Identify time step index for the start and end of the training period
option.training_start_idx=find(abs(data.timestamps-data.timestamps(1)-option.training_period(1)+1)==min(abs(data.timestamps-data.timestamps(1)-option.training_period(1)+1)),1,'first');
option.training_end_idx=find(abs(data.timestamps-data.timestamps(1)-option.training_period(2)+1)==min(abs(data.timestamps-data.timestamps(1)-option.training_period(2)+1)),1,'first');
option.nb_data_training=option.training_end_idx-option.training_start_idx;
%% Check data
if any(diff(data.timestamps)<=0)
   error(' Check your dataset: data.timestamps is not sorted in chronoligical order ') 
end


%% Construct the BDLM matrices
model=matrices_construction(model,data);

%% model option
if ~isfield(option,'method')
    option.method='kalman'; %defaut estimation method -> kalman filter
end

%% computing method
if ~isfield(option,'computing_method')
    option.computing_method='serial'; % defaut Expectation_Maximization algorithm -> serial computing
end

%% plot options
if ~isfield(option,'secondary_plots')
    option.secondary_plots=1;
end

%% Export TIKZ figure
if ~isfield(option,'export_plots')
    option.export_plots=0;
end

%% Select linewidth
if ~isfield(option,'linewidth')
    option.linewidth=1;
end

%% Select subsamples in order to reduce the number of points
if ~isfield(option,'subsample')
    option.subsample=1;
end

%% number of x-axis division
if ~isfield(option,'ndivx')
    option.ndivx=5;
end

%% number of y-axis division
if ~isfield(option,'ndivy')
    option.ndivy=3;
end
%% online with EM procedures

if ~isfield(option,'onl_estim_window')
     option.onl_estim_window=30; %day
end
if ~isfield(option,'onl_plot') % to meet the condition in BDLM.m (options 3 & 4)
    option.onl_plot=0;
end
 



