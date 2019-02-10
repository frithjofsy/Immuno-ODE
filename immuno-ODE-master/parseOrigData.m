% JN Kather 2018
% this script will pre-process temporal data from Le et al, NEJM
% source: radiographic data for CRC patients at different time points
% 
clear variables, format compact, close all, clc

rng('default');
% read the orginal data
entity = 'RES2';
spiderdata = ...
    readtable('Dataset 1(1).csv');
targetName = ['timelines_',entity];
dosave = true;
% extract data
allIDs = spiderdata.ID;
allDays = spiderdata.DAYS;
allValues = spiderdata.CHANGE_PERCENT;

% parse spider plot data 

for i=1:size(allIDs,1)
    disp(['reading row ',num2str(i),' with ID ',char(allIDs(i))]);
    if (i==1)
        disp('-> This is the first row. Starting first patient');
        patientcount = 1;
        currPatientData =  [0,0];
        collectError(patientcount,:) =  [allDays(i),allValues(i)]; 
        lastID = allIDs{i};
        
    elseif (~isempty(allIDs{i}))
        disp('-> This is the first row of the next patient. Will save old patient and start a new one in the next round.');
        patientCollection{patientcount} = currPatientData; % save data
        patientNames{patientcount} = lastID; % save patient name
        patientcount = patientcount+1; % next patient
        currPatientData = [0,0];
        collectError(patientcount,:) =  [allDays(i),allValues(i)]; 
        lastID = allIDs{i};
    elseif (i==size(allIDs,1))
        disp('-> This was the last row. Will save old patient.');
        patientCollection{patientcount} = currPatientData; % save data
        patientNames{patientcount} = lastID; % save patient name
    else
        disp('* This row is not a new patient. Concatenating data.');
        if allValues(i)<-100, allValues(i) = -100; end % less than -100% change is forbidden
        currPatientData = [currPatientData; allDays(i),allValues(i)];
    end
end

% remove all patients with less than 4 data points
rmo = [];
for i = 1:numel(patientCollection)
    if size(patientCollection{i},1)<4
       rmo = [rmo,i];
    end
end
patientCollection(rmo) = [];
patientNames(rmo) = [];

figure, hold on
allColz = lines(numel(patientNames));
% now plot spider plots
for i=1:numel(patientNames)
    currPatientName = patientNames{i};
    if contains(currPatientName,'PD')
        col = [237 28 36];
    elseif contains(currPatientName,'SD')
        col = [58 82 163];
    elseif contains(currPatientName,'PR')
        col = [0 128 0];
    else
        col = 255*allColz(i,:);
    end
    if ~isempty(patientCollection{i})
    plot(patientCollection{i}(:,1),patientCollection{i}(:,2),'+-','LineWidth',1.2,'Color',col/255)
    text(patientCollection{i}(end,1)+5,patientCollection{i}(end,2)+rand()-0.5,currPatientName,...
        'FontSize',8,'VerticalAlignment','middle');
    end
end

% decorations
axis square
set(gcf,'Color','w');

figure
boxplot(collectError,{'time (days)','size (% points)'});
title('error in data picking');

axis square
set(gcf,'Color','w');

%**************************
targetName
accuracy = round(median(collectError),2)
precision = round(iqr(collectError),2)
numPatients = numel(patientCollection)
%**************************

if dosave
save(['./Data files/',targetName,'.mat'],'patientNames','patientCollection');
end