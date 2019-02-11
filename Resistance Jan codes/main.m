clear variables;

%defining necessary settings
settings.whichModel = 1;
settings.fileName = 'dataExample.xlsx';


patients = readData(settings.fileName);
model = selectModel(settings.whichModel);

results = fitModel(patients, model);

%% 
plotResults(results);