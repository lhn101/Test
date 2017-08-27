function plot_online_estimations(estim,data,model,option,loop_input)
%% Export TIKZ figure
if ~isfield(option,'export_plots')
    option.export_plots=0;
end
%% Online components
if ~isfield(option,'onl_components')
    option.onl_components={'x^{LT}';'x^{LTc}'};
%    option.onl_components={'x^{LL}'; 'x^{LTc}'};
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


dataset_x=estim.x;
dataset_V=estim.V;
Pr_M=estim.S;
label_file='SKF';

dataset_x=real(dataset_x);
dataset_V=real(dataset_V);

plot_timeN=data.timestamps;
plot_time_1=1:option.training_end_idx;
plot_time_ex=plot_timeN(plot_time_1);
% minx=plot_time_ex(1);
% maxx=plot_time_ex(end);
minx=plot_timeN(1);
maxx=plot_timeN(end);
color_babyblue=[0 0.447 0.741];
mult_factor=4;
formatOut = 'yy-mm-dd';
date_start=datestr(plot_timeN(option.training_start_idx),formatOut);
date_prev_end=datestr(plot_timeN(option.training_end_idx),formatOut);
text_plot={date_start,date_prev_end};
Xaxis_lag=150;

idx_onl_components=find(contains(model.hidden_states_names{1}(:,1)',option.onl_components)); % Identify the hidden components to be plotted
if model.nb_class>1
    idx_prS=1;
else
    idx_prS=0;
end
nb_plots=length(idx_onl_components)+1+idx_prS;

%% Plot
close(gcf)
figure1 = figure('Color',[1 1 1], 'Position',[20 30 1250 680]);
clf

    %% State probability
if idx_prS
    subplot(nb_plots,1,idx_prS,'align')
    plot(plot_timeN(plot_time_1),Pr_M(plot_time_1,1),'color',[0.0 0.6 0.2],'Linewidth',option.linewidth*2)
    hold on
    plot(plot_timeN(plot_time_1),Pr_M(plot_time_1,2),'r','Linewidth',option.linewidth*2)
    
    % Time reference
    plot([plot_timeN(option.training_start_idx), plot_timeN(option.training_start_idx)],[0 1],...
         'color',color_babyblue,'linewidth',option.linewidth/2)
    plot([plot_timeN(option.training_end_idx), plot_timeN(option.training_end_idx)],[0 1],...
         'color',color_babyblue,'linewidth',option.linewidth/2)
    text([plot_timeN(option.training_start_idx)+5, plot_timeN(option.training_end_idx)+5],...
         [0.8 0.8],...
         text_plot,...
         'FontSize',18);
    set(gca,   'XTick'              ,linspace(plot_timeN(1),plot_timeN(end),option.ndivx)...
           ,   'XTickLabel'         ,[]...
           ,   'YTick'              ,linspace(0,1,option.ndivy)...
           ,   'Ylim'               ,[0,1]);
    ylabel('$Pr(S)$','Interpreter','Latex')
    xlim([minx-Xaxis_lag maxx])
    hold off;
end

    %% Observation -> The main response needs to be set at the first place !!!
subplot(nb_plots,1,1+idx_prS,'align')
xpl=estim.y(1,plot_time_1);
spl=estim.Vy(1,plot_time_1);
px=[plot_timeN(plot_time_1); flipud(plot_timeN(plot_time_1))]'; 
py=[xpl-sqrt(spl), fliplr(xpl+sqrt(spl))];
ypl=data.values(plot_time_1,1)';
sv = sqrt(model.R{1}(model.parameter,data.timestamps(1),data.dt_steps(1)));
psy=[ypl-sv(1,1), fliplr(ypl+sv(1,1))];

mean_xpl=nanmean(xpl(round(0.05*length(xpl)):end));
std_xpl=nanstd(xpl(round(0.1*length(xpl)):end));
miny=mean_xpl-mult_factor*std_xpl;
maxy=mean_xpl+mult_factor*std_xpl;

patch(px,psy,'red','EdgeColor','none','FaceColor','red','FaceAlpha',0.2);
hold on
plot(plot_timeN(plot_time_1),ypl','-r','Linewidth',option.linewidth/2)
plot(plot_timeN(plot_time_1),xpl,'k','Linewidth',option.linewidth)
patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

% Time references
plot([plot_timeN(option.training_start_idx), plot_timeN(option.training_start_idx)],[miny maxy],...
     'color',color_babyblue,'linewidth',option.linewidth/2)
plot([plot_timeN(option.training_end_idx), plot_timeN(option.training_end_idx)],[miny maxy],... 
     'color',color_babyblue,'linewidth',option.linewidth/2)
set(gca, 'XTick'                ,linspace(plot_timeN(1),plot_timeN(end),option.ndivx)...
       , 'XTickLabel'           ,[]...
       , 'YTick'                ,linspace(miny,maxy,option.ndivy)...
       , 'YLim'                 ,[miny,maxy]...
       );
ylabel(data.labels{1})
xlim([minx-Xaxis_lag maxx])
hold off

    %% Hidden states
loop=0;
for idx=1+idx_prS+1:nb_plots
    loop=loop+1;
    subplot(nb_plots,1,idx,'align')
    xpl=dataset_x(idx_onl_components(loop),plot_time_1);
    spl=dataset_V(idx_onl_components(loop),plot_time_1);
    px=[plot_timeN(plot_time_1); flipud(plot_timeN(plot_time_1))]'; % make closed patch
    py=[xpl-sqrt(spl) fliplr(xpl+sqrt(spl))];
    hold on
    patch(px,py,'green','EdgeColor','none','FaceColor','green','FaceAlpha',0.2);
    
    plot(plot_timeN(plot_time_1),xpl,'k','Linewidth',option.linewidth)
    
    if isfield(data,'ref') % Ref data
    plot(plot_timeN(plot_time_1),data.ref(plot_time_1,idx_onl_components(loop)),'--r')
    end
   
    mean_xpl=nanmean(xpl(round(0.05*length(xpl)):end));
    std_xpl=nanstd(xpl(round(0.05*length(xpl)):end));
    mean_spl=nanmean(sqrt(spl(round(0.05*length(xpl)):end)));
    
    miny=mean_xpl-mult_factor*(std_xpl+mean_spl);
    maxy=mean_xpl+mult_factor*(std_xpl+mean_spl);
    
    %Time reference
    plot([plot_timeN(option.training_start_idx), plot_timeN(option.training_start_idx)],[miny maxy],...
         'color',color_babyblue,'linewidth',option.linewidth/2)
    plot([plot_timeN(option.training_end_idx), plot_timeN(option.training_end_idx)],[miny maxy],...
         'color',color_babyblue,'linewidth',option.linewidth/2)    
    set(gca, 'XTick'                ,linspace(plot_timeN(1),plot_timeN(end),option.ndivx)...
           , 'XTickLabel'           ,[]...
           , 'YTick'                ,linspace(miny,maxy,option.ndivy)...
           , 'YLim'                 ,[miny,maxy]...
           );
    
    ylabel(['$' model.hidden_states_names{1}{idx_onl_components(loop),1} '$ [' '$' model.hidden_states_names{1}{idx_onl_components(loop),3}(1)...
            '$,' '$M_' model.hidden_states_names{1}{idx_onl_components(loop),2} ']$' ],'Interpreter','Latex')

end
set(gca,'XTick',linspace(plot_timeN(1),plot_timeN(end),option.ndivx));
datetick('x','yy-mm-dd','keepticks')
xlabel('Time [YY-MM-DD]')
xlim([minx-Xaxis_lag maxx])
hold off

%% Export to TIKZ
if option.export_plots
    figurename=['BDLM_' option.name '_' label_file '_no_' num2str(loop_input)];
    opts=[ 'x label style={font=\huge},',...
           'y label style={font=\huge},',...
           'legend style={font=\Large},',...
        ];
    matlab2tikz('figurehandle',figure1,'filename',['figures_OC/' figurename '.tex'] ,'standalone', true,'showInfo', false,...
                'floatFormat','%.8g','extraTikzpictureOptions','font=\huge','extraaxisoptions',opts,'width','2.4\textwidth');
end
drawnow
end