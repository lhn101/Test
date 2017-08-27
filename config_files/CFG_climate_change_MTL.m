%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file reference
option.name='Climate_change_MTL';
option.training_period=[1,60*365];% Model training start and end points [Days]

%% Data
data.labels={'Temperature'};

load('temperature_data_MTL.mat')
data.timestamps=timestamps;
data.values=temperature_data;
clear timestamps temperature_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDLM Component reference numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11: Local level
% 12: Local trend
% 13: Local acceleration
% 21: Local level compatible with local trend
% 22: Local level compatible with local acceleration
% 23: Local trend compatible with local acceleration
% 31: Perodic
% 41: Autoregressive
% 51: Dynamic regression with hidden component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model components | Choose which component to employ to build the model
model.components.block{1}={[12 31 41]};

%% Model component constrains | Take the same parameter as model class #1 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[]};