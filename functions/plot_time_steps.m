function plot_time_steps(data,option)
plot_timeN=data.timestamps;
plot_time=datevec(plot_timeN);

ndiv=8; % number of division for xlabels

figure
semilogy(plot_timeN,data.dt_steps*24)


time=round(linspace(1,size(plot_timeN,1),ndiv));
set(gca,'XTick',plot_timeN(time))
time_label=cell(1,ndiv);
for j=1:ndiv
    year=num2str(plot_time(time(j),1));
    time_label{1,j}=[year(3:4) '/' num2str(plot_time(time(j),2))];
end
set(gca,'XTickLabel',time_label)
ylabel('Time step size [h]')
xlim([min(plot_timeN),max(plot_timeN)])
xlabel('YY/MM')
