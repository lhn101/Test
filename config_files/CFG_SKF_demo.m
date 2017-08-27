%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dataset : Demo_SKF
% Author : James-A. Goulet
% Date : May. 20th 2017
% Last update : Nov. 11th 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option.name='Demo_SKF';
option.training_period=[1,2*365];% Model training start and end points [Days]

%% Data
data=load_demo_SKF_data;
data.labels={'Temp[C], y_t'};

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
model.components.block{1}={[21]};
model.components.block{2}={[12]};

%% Model component constrains | Take the same parameter as model class #1 
model.components.const{2}={[0]}; 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[]};

model.param_properties={
	'\sigma_w',	 'LcT',	 '1',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#1
	'\sigma_v',	 '',	 '1',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#2
	'\sigma_w',	 'LT',	 '2',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#3
	'\sigma_v',	 '',	 '2',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#4
'\sigma_w(11)',	'LcT',	 '12',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#5
'\sigma_w(22)',	'LcT',	 '12',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#6
'\sigma_w(11)',	'LT',	 '21',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#7
'\sigma_w(22)',	'LT',	 '21',	 'Temp[C], y_t',	 [ 0    , Inf  ]	 %#8
	'\pi_z'  ,	 '',	 '12',               '',	 [ nan  , nan  ]	 %#9
	'\pi_z'  ,	 '',	 '21',               '',	 [ nan  , nan  ]	 %#10
};
 
model.parameter=[
7.73E-07 	 %#1
0.0773   	 %#2
7.73E-08 	 %#3
0.0773   	 %#4
0.00773      %#5
0.00773      %#6
0.00773    	 %#7
0.00773    	 %#8
0.000274 	 %#9
0.000274 	 %#10
]; 
 
model.p_ref=[1   2   3   2   5   6   7   8   9  10];