option.name='SKF_Latuque';
option.training_period=[1,12*365];% Model training start and end points [Days]

%% Data
data=load_data_latuque_D(7);
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
model.components.block{1}={[23 31 31 41]};
model.components.block{2}={[13 31 31 41]};

%% Model component constrains | Take the same parameter as model class #1 
model.components.const{2}={[0 1 1 1]}; 
                     
%% Model inter-components dependence | {[model_i depends on model#]_i,[...]}
model.components.ic={[]};

model.param_properties={
	'\sigma_w',	 'TcA',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#1
	'p'       ,	 'PD1',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#2
	'\sigma_w',	 'PD1',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#3
	'p'       ,	 'PD2',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#4
	'\sigma_w',	 'PD2',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#5
	'\phi'    ,	 'AR',	 '1',	 'Displ [mm]',	 [ 0    , 1    ]	 %#6
	'\sigma_w',	 'AR',	 '1',	 'Displ [mm]',	 [ 0    , Inf  ]	 %#7
	'\sigma_v',	 '',	 '1',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#8
	'\sigma_w',	 'LA',	 '2',	 'Displ [mm]',	 [ 0    , Inf  ]	 %#9
	'\sigma_v',	 '',	 '2',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#10
	'\sigma_w(11)',	 'TcA',	 '12',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#11
	'\sigma_w(22)',	 'TcA',	 '12',	 'Displ [mm]',	 [ 0    , Inf  ]	 %#12
	'\sigma_w(33)',	 'TcA',	 '12',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#13
	'\sigma_w(11)',	 'LA',	 '21',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#14
	'\sigma_w(22)',	 'LA',	 '21',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#15
	'\sigma_w(33)',	 'LA',	 '21',	 'Displ [mm]',	 [ NaN  , NaN  ]	 %#16
	'\pi_z'   ,	 '',	 '12',	 '',	 [ NaN  , NaN  ]	 %#17
	'\pi_z'   ,	 '',	 '21',	 '',	 [ NaN  , NaN  ]	 %#18
};
 
model.parameter=[
0 	 %#1
365      	 %#2
0 	 %#3
180        	 %#4
0 	 %#5
0.75     	 %#6
0.482    	 %#7
0.3      	 %#8
4.82E-08 	 %#9
0.0482   	 %#10
0 	 %#11
4.82E-07 	 %#12
0 	 %#13
4.82E-08 	 %#14
4.82E-08 	 %#15
4.82E-08 	 %#16
0.000548 	 %#17
0.000548 	 %#18
]; 
 
model.p_ref=[1   2   3   4   5   6   7   8   9  8  11  12  13  11  12  13  17  18];