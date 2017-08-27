%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='Tamar_FLT';
option.training_period=[1,365];

%% Data
data=load_data_tamar(option.training_period);
data.labels={'Frequency [Hz]', 'Temperature [ºC]', 'Traffic load [kT]', 'Traffic pattern'};

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
model.components.block{1}={[11 41],[11 31 31 41],[11 41],[51]};

%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2 3 4],[],[4],[]};