function motion_HFcharacteristics(subinfo,TR,mot_data)
% motion_HFcharacteristics
% This script gets FD summary measures from a group of subjects and relates
% them to subject demographic caracteristics. Edit the specific
% characteristics depending on what is available in each dataset
%
% subinfo = structure, with each field containing demographic info on subs
% TR = TR in seconds
% mot_data = N sub cells, each with a time X motparam matrix
%   note1: sub order needs to be the same as in subinfo or code needs to
%       be modified to reorder
%   note2: motparam expected to be in {'x','y','z','pitch','yaw','roll'}
%       order. Code or values will need to be edited to make consistent.
%   note3: motparam rotations (and translations) expected to be in mm. Will
%       need to be converted to mm (assuming 50mm head size) if not. This
%       will affect FD calculations.
%
% Script outputs several plots and saves some data to mat files in output
% directory.
%
% CGratton

calc_FDparams = 0; % this takes a while, so turn off if running repeatedly

%%% Directory information
currDir = pwd;
outDir = sprintf('%s/HFmotion_characteristics/',currDir); % make a directory to save out info
if ~exist(outDir)
    mkdir(outDir);
end


%%%%% SUBJECT INFORMATION
% get group information for subjects
group = subinfo.Diag; %PD or CTRL
PDsubs = strcmp(group,'PD');
HCsubs = strcmp(group,'CTRL');

%%%%% MOTION INFORMATION: calculate/load HF motion data for each subject
if calc_FDparams
    FDparams = getFD(mot_data,TR,outDir);
    save([outDir 'FDparams.mat'],'FDparams','subinfo');
else
    load([outDir 'FDparams.mat']);
end


%%%% Relationship between motion and subject information

% PD vs. HC
%plot_grouphists_motion_by_characteristic(PDsubs,HCsubs,FDparams,outDir,'PD','HC');

% M vs. F
sex = subinfo.sex;
Msubs = strcmp(sex,'M'); 
Fsubs = strcmp(sex,'F');
%plot_grouphists_motion_by_characteristic(Msubs,Fsubs,FDparams,outDir,'Male','Female');

% HF pwr vs. Age
age = subinfo.ageatMR;
age(age == 0) = nan; % change "0" age to nans
%plot_correlations_motion_by_characteristic(age,FDparams,outDir,'age');

% Weight
weight = subinfo.weight;
%plot_correlations_motion_by_characteristic(weight,FDparams,outDir,'weight');

% Run ANOVA on just this subset
goodsubs = logical(~isnan(age) .* ~isnan(Msubs));
Xvars_AN = {zscore(PDsubs(goodsubs)),...
    zscore(age(goodsubs)),...
    zscore(Fsubs(goodsubs))}; 
 
Yvar = FDparams.pwr_high_rel(goodsubs,2); %y relative hf pwr

[p tbl stats terms] = anovan(Yvar,Xvars_AN,'model','interaction',...
    'varnames',{'diagnosis','age','gender'},'continuous',[2]);




end




function FDparams = getFD(mot_data_all,TR,outDir)

outDir_new = [outDir 'mot_figs/'];
if ~exist(outDir_new)
    mkdir(outDir_new);
end

make_plots = 1; % make plots
if ~make_plots
    display('Not making motion and grayplots');
end

for s = 1:length(mot_data_all)
    
    %mot_data = QC.MVM;
    mot_data = mot_data_all{s};
    
    [FDparams.pwr_high(s,:), FDparams.pwr_high_rel(s,:), FDparams.pwr_high_rel_max(s,:), FDparams.pd_full{s}, FDparams.pd_prop_full{s}, FDparams.pd_relative_full{s}, FDparams.freq{s}] = mot_FFT(mot_data,TR,make_plots);

    % save and close figure
    if make_plots
        save_fig(gcf,[outDir_new subject_list_final{s} '_mot_params.pdf']);
        close(gcf);
    end
    
    % calcualte FD (or get it from QC structure) - standard, unfiltered
    % version
    ddt_mvm_filt = [zeros(1,6); diff(mot_data)];
    FD = (sum(abs(ddt_mvm_filt),2));

    
    % filter motion params
    [FDfilt, mvm_filt] = filter_motion(TR,mot_data);
    
    % rerun FFT to make plot of filtered motion parameters
    mot_FFT(mvm_filt,TR,make_plots);
    if make_plots
        save_fig(gcf,[outDir_new subject_list_final{s} '_mot_params_filt.pdf']);
        close(gcf);
    end
    
    % make grayplots in these subjects too
    % for this you'll need a QC file with various parameters
    % Otherwise, see one of the following examples:
    %   https://github.com/cgratton/Neurohackademy_Tutorial/grayplot_Neurohackademy.py
    %   https://github.com/cgratton/DartmouthMIND_tutorial/scripts/grayplot_MIND.m
    if 0
        grayplot_HFmotion(QC,3,1,mvm_filt); % always subject 1 in this case, since separate QC per sub; 3 = after nuisance regression
        save_fig(gcf,[outDir_new 'grayplot_' subject_list_final{s} '_resid.png']);
        close(gcf);
    end
    
    % some FD params on the filtered version
    FDparams.meanFD(s) = mean(FD);
    FDparams.meanFD_filt(s) = mean(FDfilt);
    
    % some other FD params to save
    FDparams.FDover_p2(s) = sum(FD<0.2)./length(FD);
    FDparams.filtFDover_p1(s) = sum(FDfilt<0.1)./length(FDfilt);
    FDparams.filtFDover_p08(s) = sum(FDfilt<0.08)./length(FDfilt);
    FDparams.numFrames = length(FD);
    
    % save original FD and filtered FD versions to params
    FDparams.FDorig{s} = FD;
    FDparams.FDfilt{s} = FDfilt;
end

end
