%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project file_data reference
option.name='exemple_1';
option.training_period=[1,2*365];% Model training start and end points [Days]

%% Data
data.labels={'Frequency','Temperature'};
[data.timestamps,data.values]=load_data_demo; %Load data

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
model.components.block{1}={[11 41],[11 31 31 41]};

%% Model component constrains | Take the same parameter as model class #1 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[2],[]};

option.method='UD';

model.param_properties={
	'\sigma_w',	 'LL',	 '1',	 'Frequency',	 [ 0    , Inf  ]	 %#1
	'\phi'    ,	 'AR',	 '1',	 'Frequency',	 [ 0    , 1    ]	 %#2
	'\sigma_w',	 'AR',	 '1',	 'Frequency',	 [ 0    , Inf  ]	 %#3
	'\sigma_v',	 '',	 '1',	 'Frequency',	 [ NaN  , NaN  ]	 %#4
	'\sigma_w',	 'LL',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#5
	'p'       ,	 'PD1',	 '1',	 'Temperature',	 [ NaN  , NaN  ]	 %#6
	'\sigma_w',	 'PD1',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#7
	'p'       ,	 'PD2',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#8
	'\sigma_w',	 'PD2',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#9
	'\phi'    ,	 'AR',	 '1',	 'Temperature',	 [ 0    , 1    ]	 %#10
	'\sigma_w',	 'AR',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#11
	'\sigma_v',	 '',	 '1',	 'Temperature',	 [ 0    , Inf  ]	 %#12
	'\phi'    ,	 'F|T',	 '1',	 'Temperature',	 [ -Inf , Inf  ]	 %#13
};
 
model.parameter=[
2.02E-08 	 %#1
0.75     	 %#2
0.00202  	 %#3
0.00202  	 %#4
6.95E-06 	 %#5
365      	 %#6
6.95E-07 	 %#7
1        	 %#8
6.95E-07 	 %#9
0.75     	 %#10
0.695    	 %#11
0.695    	 %#12
0.01     	 %#13
]; 
 
model.p_ref=[1   2   3  12   5   6   7   8   9  10  11  12  13];
 