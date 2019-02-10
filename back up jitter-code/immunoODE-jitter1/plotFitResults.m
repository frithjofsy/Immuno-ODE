% JN Kather 2018, Moffitt Cancer Center, Tampa FL
% and J Poleszczuk 2018, NCT Heidelberg
% jakob.kather?gmail.com
% 
%% Script loads fit results and plots them
clear variables; %clearing variables from the workspace
addpath('Subroutines'); %add necessary paths

currDname = 'timelines_KNT_CEA_20'; %name of the file with data
data = load(['./Data files/',currDname,'.mat']); % load raw data for fitting
fitResults = load(['./Fit results/',currDname,'.mat']); % load raw data for fitting

%% get solutions
[~, solutions] = lossfun(fitResults.Bfinal, fitResults.patients, ...
    fitResults.initGuess); %getting solution

% rescale parameter vector if you want to look at the parameters later
[B,params] = rescaleParamVector(fitResults.Bfinal,fitResults.initGuess,numel(fitResults.patients));

% plot the fits
cnst.numRow = 4; % rows for subplot
cnst.dosave = 0; % save resulting images
cnst.dospider = 0; % draw the spider plot

[initCond,BCresponses] = plotAllPatientFits(fitResults.patients,solutions,cnst);

% if cnst.dosave % save resulting figure as PNG
%    % set(gcf,'Position',1000*[0.3560    0.0605    1.0015    0.7555]);
%    set(gcf,'Position',[ 650         209        1379         882]);
%     print(gcf,['./output/last_',currDname,'_5param.png'],'-dpng','-r400');
%     print('-painters',gcf,['./output/last_',currDname,'_5param.svg'],'-dsvg');
% end
