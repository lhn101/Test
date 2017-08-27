function interp_hc=hidden_component_WT(reg_parameters,timestamp)

t_ymd=datevec(timestamp);
t_year=t_ymd(1);
t_idx_0=[t_year 01 01 00 00 00];
t_ts_0=datenum(t_idx_0);
t_idx=timestamp+1-t_ts_0;
t_ref=reg_parameters(end-1:end);
reg_parameters=reg_parameters(1:end-2)';


reg_time=linspace(t_ref(1),t_ref(2),length(reg_parameters)+2);
reg_time=[reg_time(1)-fliplr(reg_time(2:end)-reg_time(1)) reg_time reg_time(end)+reg_time(2:end)-reg_time(1)];
reg_values=[-1 reg_parameters 1];
reg_values=[fliplr(reg_values(2:end)) reg_values fliplr(reg_values(1:end-1))];
interp_hc=interp1(reg_time,reg_values,t_idx,'spline');



