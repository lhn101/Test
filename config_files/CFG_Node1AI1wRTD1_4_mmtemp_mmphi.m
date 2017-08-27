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
data.labels={'Depl [mm]','Temp1','Temp4'};
data.values= [DNode2AI2 DNode6AI2 DNode7AI1];%Enter your dataset here; each column is a dataset; missing data -> NaN 
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
model.components.block{1}={[11 41] [11 31 31 41] [11 31 31 41]};

%% Model component constrains | Take the same parameter as model class ##1
%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}
model.components.ic={[2 3] [] []};

model.param_properties={
	'\sigma_w',	 'LL',	 '1',	 'Depl [mm]',	 [ NaN  , NaN  ]	 %#1
	'\phi'    ,	 'AR',	 '1',	 'Depl [mm]',	 [ 0  , 1]	 %#2
	'\sigma_w',	 'AR',	 '1',	 'Depl [mm]',	 [ -inf  , inf  ]	 %#3
	'\sigma_v',	 '',	 '1',	 'Depl [mm]',	 [ -inf  , inf  ]	 %#4
	'\phi'    ,	 'D|PC1',	 '1',	 'Depl [mm]',	 [ -inf  , inf  ]	 %#5
	'\phi'    ,	 'D|PC2',	 '1',	 'Depl [mm]',	 [ NaN  , NaN ]	 %#6
	'\sigma_w',	 'LL',	 '1',	 'Temp1',	 [ NaN  , NaN  ]	 %#7
	'p'       ,	 'PD1',	 '1',	 'Temp1',	 [ NaN  , NaN  ]	 %#8
	'\sigma_w',	 'PD1',	 '1',	 'Temp1',	 [ NaN  , NaN  ]	 %#9
	'p'       ,	 'PD2',	 '1',	 'Temp1',	 [ NaN  , NaN  ]	 %#10
	'\sigma_w',	 'PD2',	 '1',	 'Temp1',	 [ NaN  , NaN  ]	 %#11
	'\phi'    ,	 'AR',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#12
	'\sigma_w',	 'AR',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#13
	'\sigma_v',	 '',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#14
	'\phi'    ,	 'D|T(PD)',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#15
	'\phi'    ,	 'D|T(PD)',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#16
	'\phi'    ,	 'D|T(AR)',	 '1',	 'Temp1',	 [ NaN , NaN  ]	 %#17
	'\sigma_w',	 'LL',	 '1',	 'Temp4',	 [ NaN  , NaN  ]	 %#18
	'p'       ,	 'PD3',	 '1',	 'Temp4',	 [ NaN  , NaN  ]	 %#19
	'\sigma_w',	 'PD3',	 '1',	 'Temp4',	 [ NaN  , NaN  ]	 %#20
	'p'       ,	 'PD4',	 '1',	 'Temp4',	 [ NaN  , NaN  ]	 %#21
	'\sigma_w',	 'PD4',	 '1',	 'Temp4',	 [ NaN  , NaN  ]	 %#22
	'\phi'    ,	 'AR',	 '1',	 'Temp4',	 [ NaN , NaN    ]	 %#23
	'\sigma_w',	 'AR',	 '1',	 'Temp4',	 [ NaN , NaN  ]	 %#24
	'\sigma_v',	 '',	 '1',	 'Temp4',	 [ NaN , NaN  ]	 %#25
	'\phi'    ,	 'D|T(PD)',	 '1',	 'Temp4',	 [ NaN , NaN  ]	 %#26
	'\phi'    ,	 'D|T(PD)',	 '1',	 'Temp4',	 [ NaN , NaN  ]	 %#27
	'\phi'    ,	 'D|T(AR)',	 '1',	 'Temp4',	 [ NaN , NaN  ]	 %#28
};
model.parameter=[
0        	 %#1
0.99633  	 %#2
0.0031979 	 %#3
1.6165E-06 	 %#4
1       %#5
0      	 %#6
0        	 %#7
365.24   	 %#8
0        	 %#9
1        	 %#10
0        	 %#11
0.99618  	 %#12
0.38245  	 %#13
0.072765 	 %#14
1%-0.0016413 	 %#15
1%0.2025   	 %#16
1%0.00089471 	 %#17
0        	 %#18
365.24   	 %#19
0        	 %#20
1        	 %#21
0        	 %#22
0.99854  	 %#23
0.24287  	 %#24
0.058554 	 %#25
1%-0.0021237 	 %#26
1%0.064235 	 %#27
1%0.0006466 	 %#28
];  

 
model.p_ref=[1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18   19  20  21  22  23  24  25  26 27 28];
 %% Model PCA dependence
 PCA=[ 0.7117     -0.70249
      0.70249       0.7117];
 EYE= [1 , 0 ; 0, 1];
  model.components.PCA={PCA [] []};
% model.components.PCA={EYE [] []};

% LL =108609.7244 avec PCA
