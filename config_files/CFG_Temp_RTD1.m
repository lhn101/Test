%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='Temp_RTD1';
option.training_period=[1,2*365];

%% Data
load('P13903_Donnees.mat');  %Données de P13903_imported_data.m enregistrées
data.labels={'Temp'};
data.values= [data.Node6.AI2 ];%Enter your dataset here; each column is a dataset; missing data -> NaN 
data.timestamps= datenum(data.timestamp); %Enter the time stamps (Unix format); column vector; 


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
model.components.block{1}={[11 31 31 41]};

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[ ]};

model.param_properties={
	'\sigma_w',	 'LL',	 '1',	 'Temp',	 [ 0    , Inf  ]	 %#1
	'p'       ,	 'PD1',	 '1',	 'Temp',	 [ nan  , nan  ]	 %#2
	'\sigma_w',	 'PD1',	 '1',	 'Temp',	 [ 0    , Inf  ]	 %#3
	'p'       ,	 'PD2',	 '1',	 'Temp',	 [ nan  , nan  ]	 %#4
	'\sigma_w',	 'PD2',	 '1',	 'Temp',	 [ 0    , Inf  ]	 %#5
	'\phi'    ,	 'AR',	 '1',	 'Temp',	 [ 0    , 1    ]	 %#6
	'\sigma_w',	 'AR',	 '1',	 'Temp',	 [ 0    , Inf  ]	 %#7
	'\sigma_v',	 '',	 '1',	 'Temp',	 [ 0    , Inf  ]	 %#8
};
 
model.parameter=[
1.27E-05 	 %#1
365.2224     %#2
1.27E-06 	 %#3
1        	 %#4
1.27E-06 	 %#5
0.75     	 %#6
1.27     	 %#7
1        	 %#8
]; 
 
model.p_ref=[1  2  3  4  5  6  7  8];
