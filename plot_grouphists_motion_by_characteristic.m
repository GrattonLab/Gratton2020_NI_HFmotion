function plot_grouphists_motion_by_characteristic(var1subs,var2subs,FDparams,outDir,varname1,varname2)

% assumed order - change if not correct
motion_dirs = [1:6]
motion_directions = {'x','y','z','pitch','yaw','roll'};

for i = motion_dirs
    
    motion_dir = motion_dirs(i);
    
    % HF pwr
    tmpvals.(varname1) = FDparams.pwr_high(var1subs,motion_dir);
    tmpvals.(varname2) = FDparams.pwr_high(var2subs,motion_dir);
    [h p ci stats] = ttest2(tmpvals.(varname1),tmpvals.(varname2));
    figure; nhist(tmpvals,'smooth','proportion','samebins'); %nhist = matlab central function
    xlabel(['power in HF, ' motion_directions{motion_dir}]);
    title(sprintf('%s vs. %s: t(%d)=%.02f, p=%.03f',varname1,varname2,stats.df,stats.tstat,p))
    save_fig(gcf,[outDir varname1 'vs' varname2 '_' motion_directions{motion_dir} 'HFpwr.pdf']);
    
    % rel HF pwr
    tmpvals.(varname1) = FDparams.pwr_high_rel(var1subs,motion_dir);
    tmpvals.(varname2) = FDparams.pwr_high_rel(var2subs,motion_dir);
    [h p ci stats] = ttest2(tmpvals.(varname1),tmpvals.(varname2));
    figure; nhist(tmpvals,'smooth','proportion','samebins');
    xlabel(['relative power in HF, '  motion_directions{motion_dir}]);
    title(sprintf('%s vs. %s: t(%d)=%.02f, p=%.03f',varname1,varname2,stats.df,stats.tstat,p))
    save_fig(gcf,[outDir varname1 'vs' varname2 '_' motion_directions{motion_dir} 'relHFpwr.pdf']);
    
    % max rel HF pwr
    tmpvals.(varname1) = FDparams.pwr_high_rel_max(var1subs,motion_dir);
    tmpvals.(varname2) = FDparams.pwr_high_rel_max(var2subs,motion_dir);
    [h p ci stats] = ttest2(tmpvals.(varname1),tmpvals.(varname2));
    figure; nhist(tmpvals,'smooth','proportion','samebins');
    xlabel(['max relative power in HF, '  motion_directions{motion_dir}]);
    title(sprintf('%s vs. %s: t(%d)=%.02f, p=%.03f',varname1,varname2,stats.df,stats.tstat,p))
    save_fig(gcf,[outDir varname1 'vs' varname2 '_' motion_directions{motion_dir} 'relHFpwrMax.pdf']);
    
end

end