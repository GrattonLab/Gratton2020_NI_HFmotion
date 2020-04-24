function [pwr_high, pwr_high_rel, pwr_high_rel_max, pd, pd_prop, pd_relative, freq] = mot_FFT(mot_data,TR,varargin)
% mot_FFT(mot_data,TR,varargin)
%
% mot_data = time X motionparam (assumes 6)
% TR = TR in seconds
% varargin = info about whether or not to make extra plots (if 1 = plot,
% default is not to plot)
%
% look at FD_analysis.m for original version of this function with further
% testing of different options
%
% C Gratton

if ~isempty(varargin)
    plot_results = varargin{1};
else
    plot_results = 0;
end

Fs = 1/TR; % Sampling frequency

% PMTM approach (as used in D. Fair 2019 NI paper)
for x = 1:6
    [pd(:,x) freq] = pmtm(mot_data(:,x),8,[],Fs);

    % store pd as just a percentage
    pd_prop(:,x) = pd(:,x)./sum(pd(:,x));
    
    pd_scaled(:,x) = 10.*log10(pd(:,x));
    pd_normal(:,x) = zscore(pd_scaled(:,x));
    
    % convert to percentiles
    thisDir = pd_normal(:,x);
    pd_relative(:,x) = (thisDir - min(thisDir))./(max(thisDir) - min(thisDir))*100;
end


% Calculate HF motion
inds = freq > 0.1; % indices of freq > 0.1 Hz.
pwr_high = sum(pd(inds,:),1)./sum(pd,1); %as a proportion of total pwr
pwr_high_rel = sum(pd_relative(inds,:),1)./sum(pd_relative,1); %as a proportion of total pwr
pwr_high_rel_max = max(pd_relative(inds,:));

%from colorblind_colormap function on Matlab exchange
motion_colors = [ 1.0000         0         0; %red
    0         0    1.0000; % blue
   0.6602    0.6602    0.6602; % gray
   0.85    0.85         0; %originally yellow
        0         0         0; % black
   1.0000    0.6445         0]; % orange


if plot_results
    
    % plot results
    figure('Position',[0 0 1000 800]); 
    subplot(2,2,1:2); hold on;
    for m = 1:6
        plot(mot_data(:,m),'color',motion_colors(m,:));
    end
    xlabel('TR');
    ylabel('mm');
    title(['motion data']);
    
    % Frequency per motion parameter
    subplot(2,2,3); hold on;
    for m = 1:6
        plot(freq,pd(:,m),'color',motion_colors(m,:),'LineWidth',1.5);
    end
    ylabel('Power');
    ylim([0,0.5]);
    box off;
    xlabel('Freq (Hz)');
    title(sprintf('Power Spectrum, y HF: %.02f, Sum HF: %.02f',pwr_high(2),sum(pwr_high)));
    
    % Relative frequency
    subplot(2,2,4); hold on;
    for m = 1:6
        plot(freq,pd_relative(:,m),'color',motion_colors(m,:),'LineWidth',1.5); 
    end
    ylabel('Relative Power');
    ylim([0 100]);
    box off;
    xlabel('Freq (Hz)');
    title(sprintf('Relative Power, y HF: %.02f, All HF: %.02f, y HFmax: %.02f',pwr_high_rel(2),mean(pwr_high_rel),pwr_high_rel_max(2)));
end



end