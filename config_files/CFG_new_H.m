%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='exemple_1';
option.training_period=[1,2*365];% Model training start and end points [Days]

%% Data
data.labels={'Temperature'};
[data.timestamps,data.values]=load_data_demo_2; %Load data

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
model.components.block{1}={[11 51 41]};

%% Model component constrains | Take the same parameter as model class #1 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[]};

model.param_properties={
	'\sigma_w',	 'LL',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#1
	'\sigma_w',	 'DH',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#2
	'ycp1'    ,	 'DH',	 '1',	 'Temperature',	 [ 0    , 1    ]	 %#3
	'ycp2'    ,	 'DH',	 '1',	 'Temperature',	 [ 0    , 1    ]	 %#4
	'ycp3'    ,	 'DH',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#5
	'ycp4'    ,	 'DH',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#6
	'ycp5'    ,	 'DH',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#7
	'x_05'    ,	 'DH',	 '1',	 'Temperature',	 [ -Inf , Inf  ]	 %#8
	'd_x5'    ,	 'DH',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#9
	'\phi'    ,	 'AR',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#10
	'\sigma_w',	 'AR',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#11
	'\sigma_v',	 '',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#12
};
 
model.parameter=[
0        	 %#1
0        	 %#2
0.2      	 %#3
0.8      	 %#4
1        	 %#5
0.5      	 %#6
0.5      	 %#7
7.367E+05 	 %#8
365.22   	 %#9
0.99499     	 %#10
0.1  	 %#11
0.1  	 %#12
]; 
 
model.p_ref=[1   2   3   4   5   4   3   8   9  10  11  12];