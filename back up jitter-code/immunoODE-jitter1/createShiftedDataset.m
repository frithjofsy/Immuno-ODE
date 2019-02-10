clear variables

nRequire = -30; %days, require more than N time points

load './Data files/timelines_KNT_CEA_20.mat'
j=0;
for i = 1:numel(patientCollection)
    %if numel(patientCollection{i})>(nRequire*2)
        
        currPatientData = patientCollection{i};
     
            j = j+1;
            NewPatientCollectionTrain{j} = currPatientData + repmat([0 .5],size(currPatientData,1),1);    
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
save('./Data files/timelines_KNT_CEA_20_shift_5.mat','patientCollection','patientNames');
