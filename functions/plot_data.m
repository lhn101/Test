function plot_data(data,option)
%% Plot secondary close-up figure
if ~isfield(option,'secondary_plots')
    option.secondary_plots=1;
end

%% Export TIKZ figure
if ~isfield(option,'export_plots')
    option.export_plots=1;
end

%% Select linewidth
if ~isfield(option,'linewidth')
    option.linewidth=1;
end

%% Select subsamples in order to reduce the number of points
if ~isfield(option,'subsample')
    option.subsample=1;
end

%% number of x-axis division
if ~isfield(option,'ndivx')
    option.ndivx=5;
end

%% number of y-axis division
if ~isfield(option,'ndivy')
    option.ndivy=3;
end
ndiv_secondary=option.ndivy;
mult_factor=4;


plot_timeN=data.timestamps;
plot_time=datevec(plot_timeN);
time_fraction=0.641;
plot_time_2=round(time_fraction*length(plot_timeN)):round(time_fraction*length(plot_timeN))+(14/data.dt_ref);
plot_time_1=1:option.subsample:length(plot_timeN);
plot_label=data.labels;

nb_obser=size(data.values,2); %number of observations

for i=1:nb_obser
    xpl=data.values(:,i);
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 270])
    subplot(1,3,1:2,'align')
    
    
    max_xpl=real(max(xpl(48*7:end)));
    min_xpl=real(min(xpl(48*7:end)));
    
    plot(plot_timeN(plot_time_1),xpl(plot_time_1),'r','Linewidth',option.linewidth) 
    
    set(gca,'XTick',[linspace(plot_timeN(1),plot_timeN(end),option.ndivx)]);
    datetick('x','yy-mm','keepticks')
    ylabel(data.labels{i})
    xlabel('Time [YY-MM]')
    
    mean_xpl=nanmean(xpl(round(0.05*length(xpl)):end));
    std_xpl=nanstd(xpl(round(0.1*length(xpl)):end));
    
    miny=mean_xpl-mult_factor*std_xpl;
    maxy=mean_xpl+mult_factor*std_xpl;
    ylim([miny,maxy])
    set(gca,'YTick',[linspace(miny,maxy,option.ndivy)]);
    
     if option.secondary_plots==1
        subplot(1,3,3,'align')
        plot(plot_timeN(plot_time_2),xpl(plot_time_2),'r','Linewidth',option.linewidth)

        set(gca,'XTick',[linspace(plot_timeN(plot_time_2(1)),plot_timeN(plot_time_2(size(plot_timeN(plot_time_2),1))),ndiv_secondary)]);
        datetick('x','mm-dd','keepticks')
        year=datevec(plot_timeN(plot_time_2(1)));
        xlabel(['Time [' num2str(year(1)) '--MM-DD]'])        
        hold off
        ylim([miny,maxy])
        set(gca,'YTick',[]);        
    end

    if option.export_plots
        set(gcf,'Color',[1 1 1])
        figurename=['BDLM_' option.name '_DATA_no' num2str(i)];
        cleanfigure;
        opts=['scaled y ticks = false,',...
            'scaled x ticks = false,',...
            'y tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x label style={font=\huge},',...
            'y label style={font=\huge},',...
            'legend style={font=\Large},',...
            ['restrict y to domain=' num2str(miny-0.5*abs(miny)) ':' num2str(maxy+0.5*abs(maxy))]
            ];
        matlab2tikz('figurehandle',gcf,'filename',['saved_figures/' figurename '.tex'] ,'standalone', true,'showInfo', false,'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts);
    end
    
end

