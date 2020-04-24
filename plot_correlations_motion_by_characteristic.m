function plot_correlations_motion_by_characteristic(var,FDparams,outDir,varname,varargin)
% varargin = indexing for original matrix of subjects to use

if length(varargin) > 0
    goodsubs = logical(~isnan(var).*varargin{1});
    titleendstr = '_selectedsubs';
else
    goodsubs = ~isnan(var);
    titleendstr = '';
end

%motion_dir = 2; %motion direction to test: y = 2
motion_directions = {'x','y','z','pitch','yaw','roll'};
motion_dirs = [1:6];

for i = motion_dirs
    motion_dir = motion_dirs(i);
    % HF pwr
    plot_corr(var(goodsubs),FDparams.pwr_high(goodsubs,motion_dir),varname,'HFpwr');
    save_fig(gcf,[outDir varname '_' motion_directions{motion_dir} 'HFpwr' titleendstr '.pdf']);
    
    % relative HF pwr
    plot_corr(var(goodsubs),FDparams.pwr_high_rel(goodsubs,motion_dir),varname,'relHFpwr');
    save_fig(gcf,[outDir varname '_' motion_directions{motion_dir} 'relHFpwr' titleendstr '.pdf']);
    
    plot_corr_smooth(var(goodsubs),FDparams.pwr_high_rel(goodsubs,motion_dir),varname,'relHFpwr');
    save_fig(gcf,[outDir varname '_' motion_directions{motion_dir} 'relHFpwr' titleendstr '_smooth.pdf']);
    
    % max relative HF pwr
    plot_corr(var(goodsubs),FDparams.pwr_high_rel_max(goodsubs,motion_dir),varname,'relHFpwrMax');
    save_fig(gcf,[outDir varname '_' motion_directions{motion_dir} 'relHFpwrMax' titleendstr '.pdf']);
end

end


function plot_corr(xvar,yvar,xvarname,yvarname)

xmin = min(xvar);
xmax = max(xvar);
xrange = [xmin:(xmax-xmin)/100:xmax];

[r p] = corr(xvar,yvar,'type','Spearman');
betas = polyfit(xvar,yvar,1);
yvar_fit = polyval(betas,xrange);

figure; hold on;
scatter(xvar,yvar,20,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(xrange,yvar_fit,'k','LineWidth',2);
ylim([0.1 0.9]);
axis square;

xlabel(xvarname);
ylabel(yvarname);
title(sprintf('%svs%s: rho=%.03f, p=%.03f',xvarname,yvarname,r,p));

end

function plot_corr_smooth(xvar,yvar,xvarname,yvarname)


figure; hold on;

[f goodness output] = fit(xvar,yvar,'smoothingspline','SmoothingParam',1/length(xvar));
hL = plot(f,xvar,yvar);
set(hL(1),'color','k');
set(hL(2),'color','k','LineWidth',2);


ylim([0.1 0.9]);
axis square;

xlabel(xvarname);
ylabel(yvarname);
title(sprintf('smooth parameter = %.05f',output.p));

end