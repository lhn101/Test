function data=load_demo_SKF_data
%% Time
data.timestamps=[datenum([2017 01 01 00 00 00]):1:datenum([2017 12 31 23 59 59])]';
data.TimeData=datevec(data.timestamps);
t=data.timestamps-data.timestamps(1);

%% Temperature
m=10;
s_V=0.5;

anomaly_start=round(length(data.timestamps)/2);
baseline=[m*ones(1,anomaly_start),m+0.01*(data.timestamps(anomaly_start+1:end)'-data.timestamps(anomaly_start))];

%% Conversion to an uniform dataset
data.values=[baseline+normrnd(0,s_V,1,length(data.timestamps))]';

%plot(data.timestamp,baseline)
%hold on
%plot(data.timestamp,data.values)
%hold off