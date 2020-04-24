function [FDfilt mvm_filt] = filter_motion(TR,mvm)
% [FDfilt mvm_filt] = filter_motion(TR,mvm)
% a function that filters the motion parameters and calculates a filtered
% FD (FDfilt)
%
% Inputs:
% TR = TR in seconds
% mvm = time X 6 motion parameters
% 
% Outputs:
% FDfilt = array with filtered FD values
% mvm_filt = time X 6 motion parameters, filtered

% do the filtering
[butta buttb]=butter(1,0.1/(0.5/TR),'low');
pad = 100;
d = size(mvm);
temp_mot = cat(1, zeros(pad, d(2)), mvm, zeros(pad, d(2))); 
[temp_mot]=filtfilt(butta,buttb,double(temp_mot)); 
temp_mot = temp_mot(pad+1:end-pad, 1:d(2)); 
mvm_filt = temp_mot;

% For response analysis: testing flipped padded values
% mvm_flip = flipud(mvm);
% temp_mot2 = [mvm_flip(end-99:end,:); mvm; mvm_flip(1:100,:)];
% [temp_mot2]=filtfilt(butta,buttb,double(temp_mot2)); 
% temp_mot2 = temp_mot2(pad+1:end-pad, 1:d(2)); 
% mvm_filt2 = temp_mot2;
% figure; hold on;
% plot(mvm,'Color',[0.5 0.5 0.5],'Linewidth',3);
% plot(mvm_filt,'b','Linewidth',2);
% plot(mvm_filt2,'r','Linewidth',0.5);
% xlabel('TR');
% ylabel('mm');

% diff filtered params
ddt_mvm_filt = [zeros(1,6); diff(mvm_filt)];

% calculate FD
FDfilt=(sum(abs(ddt_mvm_filt),2));

end