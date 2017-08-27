%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='Node1AI1wRTD1_4_mmtemperarure_mmphi';
option.training_period=[1,3*365];
%tous les phi sont independants
%% Data
load('P13903_Donnees.mat','Dtimestamps');  %Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode2AI2');  %Capteur de fissure  - Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode6AI2');  %Capteur de temperature RTD5 - Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode7AI1');  %Capteur de temperature RTD6 - Données de P13903_imported_data.m enregistrées
%data.labels={'Depl [mm]','Temp1','Temp4'};
data.labels={'Depl [mm]','Temp1','Temp2','Temp3'};

%data.values= [DNode2AI2 DNode6AI2 DNode7AI1];%Enter your dataset here; each column is a dataset; missing data -> NaN 
data.values= [DNode2AI2 DNode6AI2  DNode6AI2  DNode6AI2 ];%Enter your dataset here; each column is a dataset; missing data -> NaN 
%data.values= [DNode2AI2 DNode7AI1];%Enter your dataset here; each column is a dataset; missing data -> NaN 


data.timestamps= datenum(Dtimestamps); %Enter the time stamps (Unix format); column vector; 

 option.method='UD';
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
model.components.block{1}={[11 41] [11 31 31 41]  [11 31 31 41]  [11 31 31 41] };

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2 3 4] [] [] []};