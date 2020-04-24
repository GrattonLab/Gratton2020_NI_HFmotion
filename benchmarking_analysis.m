function benchmarking_analysis(corrmat,FDparams,distmat,varargin)
% benchmarking_analysis(corrmat,FDparams,atlas,varargin)
% Script to run benchmarking calculations
% example call:
%
% Inputs
% corrmat = roi X roi X subject correlation matrix (fisher transformed)
% FDparams = structure with filtered and unfiltered FD info (see
%   motion_HFcharacteristics for how this is created)
% distmat = roi X roi matrix with distances between rois
% varargin = timecourses of unscrubbed data for random scrubbing analyses
%   tcmat{s}(regions X time) - cell structure, one cell per sub, with each
%   sub having a matrix that is regions X time
%
% CGratton
% Modified 20180614 ver: Modified to get data from Becky's latest run of PD processing with 143
% Modified 04.20.2020: Modified to remove dataset specific elements for
% publication


random_scrub = 0; % include random scrubbing comparison? Only do so for conditions with scrubbing
nperm = 100; % nperm for random scrubbing
if random_scrub
    tcmat = varargin{1};
end

%%% Directory information
currDir = pwd;
outDir = sprintf('%s/benchmarking/',currDir);
if ~exist(outDir)
    mkdir(outDir);
end

%%%% ROI INFORMATION
mask = ones(size(distmat)); % create a mask for the data
mask = logical(triu(mask,1));
distmat_lin = distmat(mask);


%%%%% MOTION analyses:

mvmt_types = {'meanFDpre','meanFDpre_filt'};

for m = 1:length(mvmt_types)
    
    disp(['******FD MEASURE: ' mvmt_types{m}]);
    mvmt_name = mvmt_types{m}; %'meanFDpre';
    mvmt = FDparams.(mvmt_name); %these are in original, not sorted, order
        
    % do correlation per edge
    [FD_fc_corr pset] = calc_FD_fc_corr(mvmt,corrmat,mask);
    
    % if this is a version with scrubbing, do comparisons to the same number of randomly scrubbed frames
    if random_scrub
        [FD_fc_corr_RAND pset_RAND] = calc_FD_fc_corr_RAND(mvmt,FDparams.numFrames,tcmat,nperm,mask);
    else
        FD_fc_corr_RAND = nan*ones([nperm,size(FD_fc_corr)]);
        pset_RAND = nan*ones([nperm,size(FD_fc_corr)]);
    end
    
    % number of significant QC-FC correlations (Bonferroni)
    persig = sum(FDR(pset)<0.05)/length(pset)*100; %FDR function in Matlab Central
    for p = 1:nperm
        if random_scrub
            persig_RAND(p) = sum(FDR(pset_RAND(p,:))<0.05)/length(pset_RAND(p,:))*100;
        else
            persig_RAND(p) = nan;
        end
    end
    disp('Percent significant connections QC-FC');
    fprintf('All: %.03f, RAND: %.03f\n',persig,mean(persig_RAND));
    
    % number of connections with corr > 0.4 to FD (not limited by sub #)
    numover04 = sum(FD_fc_corr>0.4)/length(FD_fc_corr)*100;
    for p = 1:nperm
        if random_scrub
            numover04_RAND(p) = sum(FD_fc_corr_RAND(p,:)>0.4)/length(FD_fc_corr_RAND(p,:))*100;
        else
            numover04_RAND(p) = nan;
        end
    end
    disp('Percent of QC-FC correlations > 0.4');
    fprintf('All: %.03f, RAND: %.03f\n',numover04,mean(numover04_RAND));
    
    % median QC-FC correlations
    medianr = median(abs(FD_fc_corr));
    for p = 1:nperm
        if random_scrub
            medianr_RAND(p) = median(abs(FD_fc_corr_RAND(p,:)));
        else
            medianr_RAND(p) = nan;
        end
    end
    disp('Median abs QC-FC');
    fprintf('All: %.03f,RAND: %.03f\n',medianr,mean(medianr_RAND));
    
    figure;
    bar([1 2],[mean(numover04_RAND) numover04 ],'k');
    set(gca,'Xtick',[1 2],'Xticklabel',{'RAND','all'});
    ylim([0 5]);
    ylabel('% of connections > 0.4');
    save_fig(gcf,[outDir 'mvmt_per_sign_diffs_' mvmt_name '.png']);
    
    % distance dependent plot
    disp('Dist dependent correlation, QC-FC');
    distdep = distanceplot(distmat_lin,FD_fc_corr,outDir,mvmt_name,'allsubs',1);
    for p =1:nperm
        if random_scrub
            distdep_RAND(p) = distanceplot(distmat_lin,FD_fc_corr_RAND(p,:,:,:),outDir,mvmt_name,'RAND',0);
        else
            distdep_RAND(p) = nan;
        end
    end
    fprintf('All: %.03f, RAND: %.03f\n',distdep,mean(distdep_RAND));
    
    %save out metrics of interest
    save([outDir 'CiricMetrics_' mvmt_types{m} '.mat'],'FD_fc_corr','FD_fc_corr_RAND',...
        'pset','pset_RAND','persig','persig_RAND',...
        'medianr','medianr_RAND',...
        'numover04','numover04_RAND',...
        'distdep','distdep_RAND');
    
    % close all
    close('all');
end

end


function distmat = calc_distances(atlas_params)

[roifiles x y z] = textread(atlas_params.atlas_file,'%s%f%f%f');

coords = [x y z];
d = pdist(coords);
distmat = squareform(d);

end

function [FD_fc_corr pvals] = calc_FD_fc_corr_RAND(mvmt,numFrames,tcmat,nperm,mask)

% for a given number of random permutations
for p = 1:nperm
    
    % for each subject
    for s = 1:length(tcmat)
        
        % get the timeseries
        tc = tcmat{s}'; %this should now be TIME x ROI for this subject        
        
        % randomly scrub the data to the same extent as real scrubbing
        randind = randperm(size(tc,1));
        tc_rand = tc(randind,:);
        rand_corrmat(:,:,s) = corrcoef(tc_rand(1:numFrames(s),:));        
    end
    
    % compute FD-fc correlations
    [FD_fc_corr(p,:) pvals(p,:)] = calc_FD_fc_corr(mvmt,rand_corrmat,mask);
    
    clear rand_corrmat tc_rand randind;
end
end

function [FD_fc_corr p] = calc_FD_fc_corr(meanFD,fc_matrix,mask)

% linearize fc matrix
for i = 1:size(fc_matrix,3)
    tmp = fc_matrix(:,:,i);
    fc_matrix_lin(i,:) = tmp(mask);
end

% loop over connections and calculate correlations
for fc = 1:size(fc_matrix_lin,2) 
    [r(fc) p(fc)] = corr(meanFD',fc_matrix_lin(:,fc));
end

% return z-transformed correlations
FD_fc_corr = atanh(r);

end

function r1 = distanceplot(distmat_lin,FD_fc_corr,outDir,mvmt_name,group,plotfig)

disp(['distance plots: ' group]);
r1 = corr(distmat_lin,FD_fc_corr');

if plotfig
figure;
scatter(distmat_lin,FD_fc_corr,'k.'); hold on;
[xx ind] = sort(distmat_lin);
%yy = smooth(distmat_lin,FD_fc_corr',0.1); %this takes a little while
%yyA = lowess([distmat_lin,FD_fc_corr'],0.1);
%plot(xx,yy(ind),'r-','linewidth',3);
%[slope offset] = linfit(FD_fc_corr,distmat_lin);
%yy = xx*slope + offset;
P = polyfit(distmat_lin,FD_fc_corr',1);
yy = xx*P(1) + P(2);
plot(xx,yy,'r-','linewidth',3);
hline_new(0,'k-',1);
xlabel('distance');
ylabel('FD-FC corr');
title(sprintf('r=%.03f',r1));
save_fig(gcf,[outDir 'distdep_' group '_' mvmt_name '.bmp']);
end

end

