function [timestamps,values]=load_data_demo

training_period=365;
plot_data=1;
%% Time
data.timestamp=[datenum([2016 01 01 00 00 00]):1/48:datenum([2017 12 31 23 59 59])]';
data.TimeData=datevec(data.timestamp);
t=data.timestamp-data.timestamp(1);

%% Temperature
% Local level
data.temperature.ll.mean=12;               %degree C

% Seasonal effect
data.temperature.S1.period=2*pi/(1);       %day
data.temperature.S1.phase_shift=8/24;      %day
data.temperature.S1.amplitude=4;           %degree C
data.temperature.S1.data=data.temperature.S1.amplitude*sin(data.temperature.S1.period*(t+data.temperature.S1.phase_shift));%degree K

data.temperature.S2.period=2*pi/(365);     %day
data.temperature.S2.phase_shift=8/12*365;  %day
data.temperature.S2.amplitude=9;           %degree C
data.temperature.S2.data=data.temperature.S2.amplitude*sin(data.temperature.S2.period*(t+data.temperature.S2.phase_shift));%degree K

% AR(1) -  Autoregressive model
data.temperature.AR.mean=0;              %degree C
data.temperature.AR.std=1;               %degree C
data.temperature.AR.samp_std=0.1;        %degree C
data.temperature.AR.phi=sqrt(1-(data.temperature.AR.samp_std^2/data.temperature.AR.std^2));
data.temperature.AR.data=zeros(length(t),1);

for i=1:length(t)
    if i==1
        data.temperature.AR.data(i)=normrnd(0,data.temperature.AR.std);
    else
        data.temperature.AR.data(i)=data.temperature.AR.phi*data.temperature.AR.data(i-1)+normrnd(0,data.temperature.AR.samp_std);
    end
end
data.temperature.AR.data=data.temperature.AR.data+data.temperature.AR.mean;

data.temperature.data_raw=data.temperature.ll.mean+data.temperature.S1.data+data.temperature.S2.data+data.temperature.AR.data;
data.temperature.observation_std=0.10;
data.temperature.observation = [data.temperature.data_raw+normrnd(0,data.temperature.observation_std,length(t),1)]; %DATA
data.temperature_label={'Temperature'   ;  %Sensor type
    'N/A'           ;  %Location
    'N/A'           ;  %Sensor label
    'K'            };  %Units

%% Frequencies
% Local level
data.frequency.ll.mean=1;

% Regression temperature
data.frequency.temperature_phi=0.02/std(data.temperature.data_raw);
data.frequency.temperature_reg=data.frequency.temperature_phi*(data.temperature.data_raw-data.temperature.ll.mean);

% AR(1) -  Autoregressive model
data.frequency.AR.mean=0;
data.frequency.AR.std=0.005;
data.frequency.AR.samp_std=0.0005;
data.frequency.AR.phi=sqrt(1-(data.frequency.AR.samp_std^2/data.frequency.AR.std^2));
data.frequency.AR.data=zeros(length(t),1);

for i=1:length(t)
    if i==1
        data.frequency.AR.data(i)=normrnd(0,data.frequency.AR.std);
    else
        data.frequency.AR.data(i)=data.frequency.AR.phi*data.frequency.AR.data(i-1)+normrnd(0,data.frequency.AR.samp_std);
    end
end
data.frequency.AR.data=data.frequency.AR.data+data.frequency.AR.mean;

data.frequency.data_raw=data.frequency.ll.mean+data.frequency.temperature_reg;
data.structure_response_std=0.0025; %

data.frequency.observation=[data.frequency.data_raw+normrnd(0,data.structure_response_std,length(t),1)]; %DATA
data.frequency_label={  'Frequency' ;  %Data type
    'Mode #1'   ;  %Label 1
    'N/A'       ;  %Label 2
    'Hz'       };  %Units
data.external_effects_raw=[data.temperature.data_raw];
data.external_effects_raw_labels=cat(2,data.temperature_label);
data.structure_response=[data.frequency.observation];
data.structure_response_raw=data.frequency.data_raw;
data.structure_response_labels=data.frequency_label;

%% Create a regression model for the external effects

data.DayNumber = weekday(data.timestamp);
data.TimeData=datevec(data.timestamp);
data.CumdayYear=[0 sum(eomday(2015, 1),2) sum(eomday(2015, 1:2),2) sum(eomday(2015, 1:3),2) sum(eomday(2015, 1:4),2) sum(eomday(2015, 1:5),2) sum(eomday(2015,1:6),2) sum(eomday(2015, 1:7),2) sum(eomday(2015, 1:8),2) sum(eomday(2015, 1:9),2) sum(eomday(2015, 1:10),2) sum(eomday(2015, 1:11),2) ];
data.DayYear = [data.CumdayYear(data.TimeData(:,2))'+data.TimeData(:,3)];
data.HourDay = data.TimeData(:,4)+data.TimeData(:,5)/30*0.5;
data.timestamp_idx=1:length(data.timestamp);
data.missing=any(isnan([data.HourDay data.DayYear data.DayNumber data.temperature.observation data.frequency.observation]),2);
data.not_missing_idx=find(~data.missing);
data.missing_cum=cumsum(~data.missing);
data.training_time_limit=training_period*24*2;



%% Conversion to an uniform dataset
timestamps=data.timestamp;
values=[data.frequency.observation,data.temperature.observation];

    