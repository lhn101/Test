function plot_estimations(estim,data,model,option)
%% Plot secondary close-up figure
if ~isfield(option,'secondary_plots')
    option.secondary_plots=1;
end

%% Export TIKZ figure
if ~isfield(option,'export_plots')
    option.export_plots=0;
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

dataset_x=estim.x;
dataset_V=estim.V;
Pr_M=estim.S;
label_file='SKF';

dataset_x=real(dataset_x);
dataset_V=real(dataset_V);

plot_timeN=data.timestamps;
plot_time=datevec(plot_timeN);
plot_time_init=length(plot_timeN)-(14/data.dt_ref);
plot_time_limit=length(plot_timeN);
time_fraction=0.641;
plot_time_2=round(time_fraction*length(plot_timeN)):round(time_fraction*length(plot_timeN))+(14/data.dt_ref);
plot_time_1=1:option.subsample:length(plot_timeN);
%plot_time_1=sort(randperm(length(plot_timeN),round(length(plot_timeN)/subsample)));
color_babyblue=[0 0.447 0.741];
Xaxis_lag=150;

loop=0;

%% State probabilities

if model.nb_class>1
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 270])
    plot(plot_timeN(plot_time_1),Pr_M(plot_time_1,1),'color',[0.0 0.6 0.2],'Linewidth',option.linewidth*2)
    hold on
    plot(plot_timeN(plot_time_1),Pr_M(plot_time_1,2),'r','Linewidth',option.linewidth*2)
    set(gca,'XTick'            ,linspace(plot_timeN(plot_time_1(1)),plot_timeN(plot_time_1(size(plot_timeN(plot_time_1),1))),option.ndivx),...
            'YTick'            ,[0 0.5 1],...
            'box'              ,'off');
    datetick('x','yy-mm','keepticks')
    xlabel('Time [YY-MM]')
    ylabel('$Pr(S)$','Interpreter','Latex')
    legend('Normal','Abnormal')
    xlim([plot_timeN(1)-Xaxis_lag,plot_timeN(end)])
    hold off;
    if option.export_plots
        set(gcf,'Color',[1 1 1])
        figurename=['BDLM_SKF_plot_no' num2str(loop)];
        cleanfigure('targetResolution',300);
        opts_1=['scaled y ticks = false,',...
            'scaled x ticks = false,',...
            'x label style={font=\huge},',...
            'y label style={font=\huge},',...
            'legend style={font=\Large},',...
            ];
        matlab2tikz('figurehandle',gcf,'filename',['saved_figures/' figurename '.tex'] ,'standalone', true,'showInfo', false,...
                    'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts_1,'width','1.8\textwidth');
    end
    
end

if option.secondary_plots==1
    idx_supp_plot=0;
else
    idx_supp_plot=1;
end

%% Hidden states

for idx=1:size(dataset_x,1)
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 270])
    subplot(1,3,1:2+idx_supp_plot,'align')
    xpl=dataset_x(idx,plot_time_1);
    spl=dataset_V(idx,plot_time_1);
    
    mean_xpl=nanmean(xpl(round(0.05*length(xpl)):end));
    std_xpl=nanstd(xpl(round(0.05*length(xpl)):end));
    mean_spl=nanmean(sqrt(spl(round(0.05*length(xpl)):end)));
    mult_factor=6;
    
    miny=mean_xpl-mult_factor*(std_xpl+mean_spl);
    maxy=mean_xpl+mult_factor*(std_xpl+mean_spl);
    px=[plot_timeN(plot_time_1); flipud(plot_timeN(plot_time_1))]'; % make closed patch
    py=[xpl-sqrt(spl) fliplr(xpl+sqrt(spl))];
    hold on
    patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);
    plot(plot_timeN(plot_time_1),xpl,'k','Linewidth',option.linewidth)
    
    if isfield(data,'ref') % Reference data
        plot(plot_timeN(plot_time_1),data.ref(plot_time_1,idx),'--r')
    end
    
    ylabel(['$' model.hidden_states_names{1}{idx,1} '$ [' '$' model.hidden_states_names{1}{idx,3}(1) '$,' '$M_' model.hidden_states_names{1}{idx,2} ']$' ],'Interpreter','Latex')
    set(gca,'XTick'                ,linspace(plot_timeN(plot_time_1(1)),plot_timeN(plot_time_1(size(plot_timeN(plot_time_1),1))),option.ndivx),...
            'YTick'                ,linspace(miny,maxy,option.ndivy),...
            'Ylim'                 ,[miny,maxy],...
            'box'                  ,'off');
    datetick('x','yy-mm','keepticks')
    xlabel('Time [YY-MM]')
    xlim([plot_timeN(1)-Xaxis_lag,plot_timeN(end)])
    hold off
  
    if loop==0
        if estim.smooth==0
            h=legend('$\mu_{t|t}\pm\sigma_{t|t}$','$\mu_{t|t}$');%,'Location','northoutside','Orientation','horizontal');
        elseif estim.smooth==1
            h=legend('$\mu_{t|T}\pm\sigma_{t|T}$','$\mu_{t|T}$');%,'Location','northoutside','Orientation','horizontal');
        end
        set(h,'Interpreter','Latex')
        PatchInLegend = findobj(h, 'type', 'patch');
        set(PatchInLegend, 'facea', 0.5)
    end
    
    if option.secondary_plots==1
        subplot(1,3,3,'align')
        xpl=estim.x(idx,plot_time_2);
        spl=estim.V(idx,plot_time_2);
        
        px=[plot_timeN(plot_time_2); flipud(plot_timeN(plot_time_2))]'; % make closed patch
        py=[xpl-sqrt(spl), fliplr(xpl+sqrt(spl))];
        hold on
        patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);
        plot(plot_timeN(plot_time_2),xpl,'k','Linewidth',option.linewidth)
        set(gca,'XTick'           ,linspace(plot_timeN(plot_time_2(1)),plot_timeN(plot_time_2(size(plot_timeN(plot_time_2),1))),ndiv_secondary),...
                'YTick'           ,[],...
                'Ylim'            ,[miny,maxy],...
                'box'             ,'off');
        datetick('x','mm-dd','keepticks')
        year=datevec(plot_timeN(plot_time_2(1)));
        xlabel(['Time [' num2str(year(1)) '--MM-DD]'])
        hold off     
    end
    
    loop=loop+1;
    if option.export_plots
        set(gcf,'Color',[1 1 1])
        figurename=['BDLM_' option.name '_' label_file '_no' num2str(loop)];
        cleanfigure('targetResolution',300);
        opts=['scaled y ticks = false,',...
            'scaled x ticks = false,',...
            'y tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x label style={font=\huge},',...
            'y label style={font=\huge},',...
            'legend style={font=\Large},',...
            ['restrict y to domain=' num2str(miny-0.5*abs(miny)) ':' num2str(maxy+0.5*abs(maxy))]
            ];
        matlab2tikz('figurehandle',gcf,'filename',['saved_figures/' figurename '.tex'] ,'standalone', true,'showInfo', false,...
                    'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts,'width','1.8\textwidth');
    end  
end

%% Observation & estimated values

for i=1:model.nb_obs
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 270])
    subplot(1,3,1:2+idx_supp_plot,'align')
    xpl=estim.y(i,plot_time_1);
    spl=estim.Vy(i,plot_time_1);
    px=[plot_timeN(plot_time_1); flipud(plot_timeN(plot_time_1))]'; % make closed patch
    py=[xpl-sqrt(spl), fliplr(xpl+sqrt(spl))];
    ypl=data.values(plot_time_1,i)';
    sv = sqrt(model.R{1}(model.parameter,data.timestamps(1),data.dt_steps(1)));
    psy=[ypl-sv(i,i), fliplr(ypl+sv(i,i))];
    
    mean_xpl=mean(xpl(round(0.25*length(xpl)):end));
    std_xpl=std(xpl(round(0.25*length(xpl)):end));
    mult_factor=6;
    miny=mean_xpl-mult_factor*std_xpl;
    maxy=mean_xpl+mult_factor*std_xpl;
    
    patch(px,psy,'red','EdgeColor','none','FaceColor','red','FaceAlpha',0.2);
    hold on
    plot(plot_timeN(plot_time_1),data.values(plot_time_1,i),'-r','Linewidth',option.linewidth) 
    patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);
    plot(plot_timeN(plot_time_1),xpl,'k','Linewidth',option.linewidth)  
    set(gca,'XTick'            ,linspace(plot_timeN(plot_time_1(1)),plot_timeN(plot_time_1(size(plot_timeN(plot_time_1),1))),option.ndivx),...
            'YTick'            ,linspace(miny,maxy,option.ndivy),...
            'Ylim'             ,[miny,maxy],...
            'box'              ,'off');
    datetick('x','yy-mm','keepticks')
    xlabel('Time [YY-MM]')
    ylabel(data.labels{i})
    xlim([plot_timeN(1)-Xaxis_lag,plot_timeN(end)])
    hold off
    
    if i==1
        if estim.smooth==0
            h=legend('$y_{t}\pm \sigma_V$','$y_{t}$','$E[Y_t|y_{1:t}]\pm\sigma_E[Y_t|y_{1:t}]$','$E[Y_t|y_{1:t}]$');%,'Location','northoutside','Orientation','horizontal');
        elseif estim.smooth==1
            h=legend('$y_{t}\pm \sigma_V$','$y_{t}$','$E[Y_t|y_{1:T}]\pm\sigma_E[Y_t|y_{1:T}]$','$E[Y_t|y_{1:T}]$');%,'Location','northoutside','Orientation','horizontal');
        end
        set(h,'Interpreter','Latex')
        PatchInLegend = findobj(h, 'type', 'patch');
        set(PatchInLegend, 'facea', 0.5)
        %legend boxoff
    end
    
    if option.secondary_plots==1
        subplot(1,3,3,'align')
        xpl=estim.y(i,plot_time_2);
        spl=estim.Vy(i,plot_time_2);
        
        px=[plot_timeN(plot_time_2); flipud(plot_timeN(plot_time_2))]'; % make closed patch
        py=[xpl-sqrt(spl), fliplr(xpl+sqrt(spl))];
        psy=[xpl-sv(i,i), fliplr(xpl+sv(i,i))];
        patch(px,psy,'red','EdgeColor','none','FaceColor','red','FaceAlpha',0.2);
        hold on
        plot(plot_timeN(plot_time_2),data.values(plot_time_2,i),'-r')
        patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);
        plot(plot_timeN(plot_time_2),xpl,'k','Linewidth',option.linewidth)
 
        set(gca,'XTick'          ,linspace(plot_timeN(plot_time_2(1)),plot_timeN(plot_time_2(size(plot_timeN(plot_time_2),1))),ndiv_secondary),...
                'YTick'          ,[],...
                'Ylim'           ,[miny,maxy],...
                'box'            ,'off');
        datetick('x','mm-dd','keepticks')
        year=datevec(plot_timeN(plot_time_2(1)));
        xlabel(['Time [' num2str(year(1)) '--MM-DD]'])
        %ylabel(option.meas_labels{i})
        hold off       
    end
    
    loop=loop+1;   
    if option.export_plots
        set(gcf,'Color',[1 1 1])           
        figurename=['BDLM_' option.name '_' label_file '_no' num2str(loop)];
        cleanfigure('targetResolution',300);
        opts=['scaled y ticks = false,',...
            'scaled x ticks = false,',...
            'y tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
            'x label style={font=\huge},',...
            'y label style={font=\huge},',...
            'legend style={font=\Large},',...
            ['restrict y to domain=' num2str(miny-0.5*abs(miny)) ':' num2str(maxy+0.5*abs(maxy))]
            ];
        matlab2tikz('figurehandle',gcf,'filename',['saved_figures/' figurename '.tex'] ,'standalone', true,'showInfo', false,...
                    'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts,'width','1.8\textwidth');
    end 
end
%% Parameters estimations

if isfield(estim,'parameter')
    bound_param=[model.param_properties{:,5}];
    bound_param=reshape(bound_param,[],length(model.p_ref))';
    idx_param_estim=find(and(~isnan(bound_param(:,1)),~isnan(bound_param(:,1))));
    for k=1:length(idx_param_estim)
        maxy=max(estim.parameter(:,idx_param_estim(k)));
        miny=min(estim.parameter(:,idx_param_estim(k)));
        FigHandle = figure;
        plot(plot_timeN(plot_time_1),estim.parameter(:,idx_param_estim(k)),'color',color_babyblue,'linewidth',2*option.linewidth)
        hold on
        set(FigHandle, 'Position', [100, 100, 1300, 270])
        set(gca,'XTick'            ,linspace(plot_timeN(plot_time_1(1)),plot_timeN(plot_time_1(size(plot_timeN(plot_time_1),1))),option.ndivx),...
                'YTick'            ,linspace(miny,maxy,option.ndivy),...
                'box'              ,'off');
        datetick('x','yy-mm','keepticks')    
        xlabel('Time [YY-MM]')
        ylabel(['$' model.param_properties{idx_param_estim(k),1}  '^{' model.param_properties{idx_param_estim(k),2} '}$'],'Interpreter','Latex')
        xlim([plot_timeN(1)-Xaxis_lag,plot_timeN(end)])
        hold off
        
        loop=loop+1;
        if option.export_plots
            set(gcf,'Color',[1 1 1])           
            figurename=['BDLM_' option.name '_' label_file '_param_no' num2str(loop)];
            cleanfigure('targetResolution',300);
            opts=['scaled y ticks = false,',...
                'scaled x ticks = false,',...
                'y tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
                'x tick label style={/pgf/number format/.cd, fixed, fixed zerofill,precision=2},',...
                'x label style={font=\huge},',...
                'y label style={font=\huge},',...
                'legend style={font=\Large},'
%                 ['restrict y to domain=' num2str(miny-0.5*abs(miny)) ':' num2str(maxy+0.5*abs(maxy))]
                ];
            matlab2tikz('figurehandle',gcf,'filename',['saved_figures/' figurename '.tex'] ,'standalone', true,'showInfo', false,...
                        'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts,'width','1.8\textwidth');
        end
    end
end