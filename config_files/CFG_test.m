%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='test';
option.training_period=[1,730];

%% Data
data.labels={'depl','temp'};
data.values= %Enter your dataset here; each column is a dataset; missing data -> NaN 
data.timestamps= %Enter the time stamps (Unix format); column vector; 

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
model.components.block{1}={[11 41 ] [11 31 31 41 ] };

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2 ] [ ] };

