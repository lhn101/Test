function plot_missing_data(data,option)
plot_timeN=data.timestamps;
plot_time=datevec(plot_timeN);

nb_obser=size(data.values,2); %number of observations
ndiv=8; % number of division for xlabels

if ~any(any(isnan(data.values)))
   disp(' ') 
      disp('    -> No data is missing.') 
   disp(' ') 
else
figure
for i=1:nb_obser
    idx=find(isnan(data.values(:,i)));
    scatter(plot_timeN(idx),idx*0+i,'xr')
    hold on
end
hold off

time=round(linspace(1,size(plot_timeN,1),ndiv));
set(gca,'XTick',plot_timeN(time))
time_label=cell(1,ndiv);
for j=1:ndiv
    year=num2str(plot_time(time(j),1));
    time_label{1,j}=[year(3:4) '/' num2str(plot_time(time(j),2))];
end
set(gca,'XTickLabel',time_label)
set(gca,'YTick',[1 2 3])
data_label=cell(1,nb_obser);
for j=1:nb_obser
    data_label{1,j}=[data.labels{j}];
end
set(gca,'yTickLabel',data_label)
xlim([min(plot_timeN),max(plot_timeN)])
xlabel('YY/MM')
end
