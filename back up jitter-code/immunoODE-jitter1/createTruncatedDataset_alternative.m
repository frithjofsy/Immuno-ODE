clear variables

%nRequire = 120; %days, require more than N time points

load './Data files/timelines_KNT_CEA_20.mat'
j=0;
for i = 1:numel(patientCollection)
    %if numel(patientCollection{i})>(nRequire*2)
        
        currPatientData = patientCollection{i};
        %NewPatientCollectionFull{j} = patientCollection{i};
        %nRequire = max(patientCollection{i}(:,1))/2;
        numT = numel(patientCollection{i}(:,1));
        indx = 1:ceil(numT*(2/3));
        %if sum(indx) >= 6
            j = j+1;
            NewPatientCollectionTrain{j} = currPatientData(indx,:);    
            %NewPatientCollectionTest{j} = currPatientData((nRequire+1):end,:);    
            NewPatientNames{j} = patientNames{i};
        %end
    %end
end
    
% prepare training cohort
clear patientCollection patientNames
patientCollection = NewPatientCollectionTrain;
patientNames = NewPatientNames;

% save training cohort
save('./Data files/timelines_KNT_Crop33perc.mat','patientCollection','patientNames');
