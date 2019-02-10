% JN Kather 2018, Moffitt Cancer Center, Tampa FL
% and J Poleszczuk 2018, NCT Heidelberg
% jakob.kather?gmail.com
% 
%% Script loads fit results for <90 days fitting, but 
%% plots the solutions until the last available data point
clear variables;
addpath(genpath('./subroutines')); %add necessary paths

currDname = 'timelines_NCT_CEA_20'; %name of the file with data
currDname2 = 'timelines_NCT_CEA_20_jitter'; 
% timelines_KNT_CEA_20_jitter OR timelines_KNT_crop100days_jitter OR timelines_NCT_CEA_20_jitter
cnst.plotVerticalLine = 0;%100; 'Crop33perc'
cnst.numRow = 4; % rows for subplot
cnst.dosave = 0; % save resulting images
cnst.loQuantile = 0.05;
data = load(['./Data files/',currDname,'.mat']); % load raw data for fitting
fitRes = load(['./Fit results/',currDname2,'.mat']); % load raw data for fitting

for jit = 1:numel(fitRes.allQs) % for each jitter experiment
disp(['****',newline,'starting jitter ',num2str(jit)]);
%% get solutions. First: reformat data
patients = struct();
for i = 1:numel(data.patientCollection) %for each patient
    patients(i).timelineBIO{1} = data.patientCollection{i};
    patients(i).newID = strrep(strrep(strrep(data.patientNames{i},'LE_CRC_',''),'_BIOCHEM',''),'_','-');
end

[~, solutions] = lossfun(fitRes.allBs(:,jit), patients, ...
    fitRes.allIs(:,jit)); % getting solution

[t1{jit},yd{jit},xi{jit},Ti{jit},Ei{jit},patnames{jit},R2{jit}] =  ...
                collectJitterFits(patients,solutions); % collect jittered fits
            
end

% now plot
plotJitterFits(t1,yd,xi,Ti,Ei,patnames,R2,fitRes.allQs,cnst); % plot the fits
