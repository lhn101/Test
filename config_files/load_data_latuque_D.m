function data_sim=load_data_latuque_D(inp)

load('Latuque_raw_data_D_T1H.mat')  %Load raw data from Latuque Dam
%% Displacement
data.timestamp={raw_data.timestamp.LTU017PIAPRG300 raw_data.timestamp.LTU017PIAPRG300 raw_data.timestamp.LTU017PIAPRG300 raw_data.timestamp.LTU016PIAPE910 raw_data.timestamp.LTU016PIAPE910 raw_data.timestamp.LTU016PIAPE910 raw_data.timestamp.LTU014PIAEVA920 raw_data.timestamp.LTU014PIAEVA920 raw_data.timestamp.LTU014PIAEVA920 raw_data.timestamp.LTU012PIAPE010 raw_data.timestamp.LTU012PIAPE010 raw_data.timestamp.LTU012PIAPE010 raw_data.timestamp.LTU010PIAPRG910 raw_data.timestamp.LTU010PIAPRG910 raw_data.timestamp.LTU010PIAPRG910};

data.displacement={raw_data.LTU017PIAPRG300.X,raw_data.LTU017PIAPRG300.Y,raw_data.LTU017PIAPRG300.Z,raw_data.LTU016PIAPE910.X,raw_data.LTU016PIAPE910.Y,raw_data.LTU016PIAPE910.Z,raw_data.LTU014PIAEVA920.X,raw_data.LTU014PIAEVA920.Y,raw_data.LTU014PIAEVA920.Z,raw_data.LTU012PIAPE010.X,raw_data.LTU012PIAPE010.Y,raw_data.LTU012PIAPE010.Z,raw_data.LTU010PIAPRG910.X,raw_data.LTU010PIAPRG910.Y,raw_data.LTU010PIAPRG910.Z};
data.displacement_label={'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement' 'Displacement';
    'LTU017PIAPRG300_X' 'LTU017PIAPRG300_Y' 'LTU017PIAPRG300_Z' 'LTU016PIAPE910_X' 'LTU016PIAPE910_Y' 'LTU016PIAPE910_Z' 'LTU014PIAEVA920_X' 'LTU014PIAEVA920_Y' 'LTU014PIAEVA920_Z' 'LTU012PIAPE010_X' 'LTU012PIAPE010_Y' 'LTU012PIAPE010_Z' 'LTU010PIAPRG910_X' 'LTU010PIAPRG910_Y' 'LTU010PIAPRG910_Z';
    'mm'                'mm'                'mm'                'mm'               'mm'               'mm'               'mm'                'mm'                'mm'                'mm'               'mm'               'mm'               'mm'                'mm'                 'mm'}; 
data.displacement_std=0.005*ones(1,1);

%% Temperature
data.timestamp_temp={raw_data.timestamp.LTU038THAPRG900 raw_data.timestamp.LTU037THAPRG900 raw_data.timestamp.LTU036THAPRG900 raw_data.timestamp.LTU035THAPRG900};
data.temperature={raw_data.LTU038THAPRG900 raw_data.LTU037THAPRG900 raw_data.LTU036THAPRG900 raw_data.LTU035THAPRG900};
% data.temperature{data.temperature>50}=NaN;
% data.temperature(data.temperature<-50)=NaN;
% data.temperature=data.temperature; %Kelvin conversion
data.temperature_std=0.04*ones(1,9);
data.temperature_label={'Temperature'     'Temperature'      'Temperature'       'Temperature';
                        'LTU038THAPRG900' 'LTU037THAPRG900'  'LTU036THAPRG900'   'LTU035THAPRG900';
                        'º C'             'º C'              'º C'               'º C'};

% disp(' ')
% disp('-------------------------------------------------------------------------------')
% disp('/ Select Displacement data')
% disp('-------------------------------------------------------------------------------')
% disp(' ')
% 
% for idx=1:length(data.displacement_label)
%     disp([' ' num2str([idx],'%02.0f\n') ' -> ' repmat([' '],1,16-length(data.displacement_label{1,idx})) data.displacement_label{1,idx} repmat([' '],1,16-length(data.displacement_label{1,idx})) data.displacement_label{2,idx} repmat([' '],1,20-length(data.displacement_label{2,idx})),repmat([' '],1,11-length(data.displacement_label{3,idx})) '[' data.displacement_label{3,idx} ']' repmat([' '],1,11-length(data.displacement_label{3,idx}))])
% end 
% inputs.inp_3 = input(['Selection [1-' num2str(length(data.displacement_label)) '] : ']);
inputs.inp_3=inp;
data.displacement=data.displacement{inputs.inp_3};
data.timestamp=data.timestamp{inputs.inp_3};

%% Conversion to an uniform dataset
data_sim.timestamps=data.timestamp;
data_sim.values=[data.displacement];