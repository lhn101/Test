option.name='Latuque';
option.training_period=[1,12*365];% Model training start and end points [Days]

%% Data
data=load_data_latuque_D(1);
data.labels={'Displ [mm]'};

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
model.components.block{1}={[13 31 41]};

%% Model component constrains | Take the same parameter as model class #1 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[]};

 