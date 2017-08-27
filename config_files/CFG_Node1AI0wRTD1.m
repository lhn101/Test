%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='Node1AI0wRTD1';
option.training_period=[1,730];

%% Data
load('P13903_Donnees.mat','Dtimestamps');  %Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode1AI0');  %Données de P13903_imported_data.m enregistrées
load('P13903_Donnees.mat','DNode6AI2');  %Données de P13903_imported_data.m enregistrées
data.values= [DNode1AI0 DNode6AI2 ];%Enter your dataset here; each column is a dataset; missing data -> NaN 
data.timestamps= [Dtimestamps]; %Enter the time stamps (Unix format); column vector 
data.labels={'Depl [mm]','Temp'};

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
model.components.block{1}={[11 41 ] [11 31 31 41] };

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2] [ ] };

model.param_properties={
	'\sigma_w',	 'LL',	 '1',	 'Depl [mm]',	 [ NaN  , NaN  ]	 %#1
	'\phi'    ,	 'AR',	 '1',	 'Depl [mm]',	 [ 0    , 1    ]	 %#2
	'\sigma_w',	 'AR',	 '1',	 'Depl [mm]',	 [ 0    , Inf  ]	 %#3
	'\sigma_v',	 '',	 '1',	 'Depl [mm]',	 [ 0    , Inf  ]	 %#4
	'\sigma_w',	 'LL',	 '1',	 'Temp',         [ NaN  , NaN  ]	 %#5  
	'p'       ,	 'PD1',	 '1',	 'Temp',    	 [ NaN  , NaN  ]	 %#6
	'\sigma_w',	 'PD1',	 '1',	 'Temp',         [ NaN  , NaN  ]	 %#7
	'p'       ,	 'PD2',	 '1',	 'Temp',         [ NaN  , NaN  ]	 %#8
	'\sigma_w',	 'PD2',	 '1',	 'Temp',         [ NaN  , NaN  ]	 %#9
	'\phi'    ,	 'AR',	 '1',	 'Temp',         [ nan  , nan  ]	 %#10        OPTIMISE DANS CFG_TEMP_RTD1
	'\sigma_w',	 'AR',	 '1',	 'Temp',         [ nan  , nan  ]	 %#11        OPTIMISE DANS CFG_TEMP_RTD1
	'\sigma_v',	 '',	 '1',	 'Temp',         [ nan  , nan  ]	 %#12        OPTIMISE DANS CFG_TEMP_RTD1
	'\phi'    ,	 'D|T',	 '1',	 'Temp',         [ -Inf , Inf   ]	 %#13
};


model.parameter=[
0        	 %#1
0.98138  	 %#2
0.024549 	 %#3
1e-3 	 %#4
0        	 %#5
365.24   	 %#6
0        	 %#7
1        	 %#8
0        	 %#9
0.99611  	 %#10
0.38612  	 %#11
0.063138 	 %#12
-0.012   	 %#13
]; 

model.p_ref=[1   2   3   4   5   6   7   8   9  10  11  12  13];