clear variables

nRequire = -30; %days, require more than N time points

load './Data files/timelines_KNT_CEA_20.mat'
j=0;
for i = 1:numel(patientCollection)
    %if numel(patientCollection{i})>(nRequire*2)
        
        currPatientData = patientCollection{i};
        %NewPatientCollectionFull{j} = patientCollection{i};
        if nRequire<0
            indx = patientCollection{i}(:,1) <= (max(patientCollection{i}(:,1))+nRequire);
        else
            indx = patientCollection{i}(:,1) <= nRequire;
        end
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
save('./Data files/timelines_KNT_cropLast30d.mat','patientCollection','patientNames');
