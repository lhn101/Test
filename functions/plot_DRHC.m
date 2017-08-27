function plot_DRHC(model,option)
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

p_idx=[];
p_idx_dep=[];

for i=1:size(model.param_properties,1)
    if strcmp(model.param_properties{i,2},'DH')
        p_idx=[p_idx,i];
        
        if isnan(model.param_properties{i,5}(1))
            p_idx_dep=[p_idx_dep,1];
        else
            p_idx_dep=[p_idx_dep,0];   
        end
    end
    
end

p_DRHC=model.parameter(model.p_ref);
p_DRHC=p_DRHC(p_idx);
p_DRHC=p_DRHC(2:end);

x=linspace(p_DRHC(end-1),p_DRHC(end-1)+p_DRHC(end),500);
interp_h=zeros(1,length(x));
for i=1:length(x)
    interp_h(i)=hidden_component(p_DRHC,x(i));
end
interp_h=interp_h-interp_h(1);



xpl=x;
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 800, 370])

px=p_DRHC(1:end-2);
plot(x,interp_h,'b','Linewidth',option.linewidth)
hold on
p_idx_dep=p_idx_dep(2:end-2);
x_sc=linspace(x(1),x(end),length(px)+2)';
x_sc=x_sc(2:end-1);
plot(x_sc(find(p_idx_dep)),px(find(p_idx_dep)),'xr','Linewidth',option.linewidth,'MarkerSize',10);
plot(x_sc(find(~p_idx_dep)),px(find(~p_idx_dep))','or','Linewidth',option.linewidth);
plot(x(1),0,'md','Linewidth',option.linewidth)
plot(x(end),0,'md','Linewidth',option.linewidth)

hold off
set(gca,'XTick',[linspace(x(1),x(end)+0.001,option.ndivx)]);
datetick('x','yy-mm','keepticks')
ylabel('Hidden covariate')
xlabel('Time [YY-MM]')

ylim([0,1])
set(gca,'YTick',[0 0.25 0.5 0.75 1]);

if option.export_plots
    set(gcf,'Color',[1 1 1])
    figurename=['BDLM_' option.name '_DRHC_no' num2str(i)];
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

