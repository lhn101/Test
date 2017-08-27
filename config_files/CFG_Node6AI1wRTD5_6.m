%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='Node6AI1wRTD5_6';
option.training_period=[1,2*365];

%% Data
load('P13903_Donnees.mat','Dtimestamps');  %Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode6AI1');  %Capteur de fissure  - Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode7AI2');  %Capteur de temperature RTD5 - Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode7AI3');  %Capteur de temperature RTD6 - Données de P13903_imported_data.m enregistrées
data.labels={'Depl [mm]','Temp5','Temp6'};
data.values= [DNode6AI1 DNode7AI2 DNode7AI3];%Enter your dataset here; each column is a dataset; missing data -> NaN 
data.timestamps= datenum(Dtimestamps); %Enter the time stamps (Unix format); column vector; 


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
model.components.block{1}={[11 41] [11 31 31 41] [11 31 31 41]};

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2 3] [] []};

%% Model PCA dependence
model.components.PCA={eye(2) [] []};


