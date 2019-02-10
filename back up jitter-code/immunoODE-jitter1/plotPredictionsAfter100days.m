% JN Kather 2018, Moffitt Cancer Center, Tampa FL
% and J Poleszczuk 2018, NCT Heidelberg
% jakob.kather?gmail.com
% 
%% Script loads fit results for <100 days fitting, but 
%% plots the solutions until the last available data point
clear variables;

addpath(genpath('./subroutines')); %add necessary paths

currDname = 'timelines_KNT_CEA_20'; %name of the file with data
currDname2 = 'timelines_KNT_crop100days';
cnst.plotVerticalLine = 100;
data = load(['./Data files/',currDname,'.mat']); % load raw data for fitting
fitResults = load(['./Fit results/',currDname2,'.mat']); % load raw data for fitting

%% get solutions
%reformat data
patients = struct();
for i = 1:numel(data.patientCollection) %for each patient
    patients(i).timelineBIO{1} = data.patientCollection{i};
    patients(i).newID = strrep(strrep(strrep(data.patientNames{i},'LE_CRC_',''),'_BIOCHEM',''),'_','-');
end

[~, solutions] = lossfun(fitResults.Bfinal, patients, ...
    fitResults.initGuess); %getting solution

% rescale parameter vector if you want to look at the parameters later
[B,params] = rescaleParamVector(fitResults.Bfinal,fitResults.initGuess,numel(fitResults.patients));

% plot the fits
cnst.numRow = 4; % rows for subplot
cnst.dosave = 0; % save resulting images
cnst.dospider = 0; % draw the spider plot

[initCond,BCresponses] = plotAllPatientFits(patients,solutions,cnst);

%    gof(i,:) = rsquared(y1,ymodelERR(:,1));
%    names{i} = res(i).newID;
%    numInterpDataPoints(i) = numel(y1);%+numel(EFF);

% if cnst.dosave % save resulting figure as PNG
%    % set(gcf,'Position',1000*[0.3560    0.0605    1.0015    0.7555]);
%    set(gcf,'Position',[ 650         209        1379         882]);
%     print(gcf,['./output/last_',currDname,'_5param.png'],'-dpng','-r400');
%     print('-painters',gcf,['./output/last_',currDname,'_5param.svg'],'-dsvg');
% end


