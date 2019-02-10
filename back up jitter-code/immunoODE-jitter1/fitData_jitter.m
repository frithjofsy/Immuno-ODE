% JN Kather 2018, Moffitt Cancer Center, Tampa FL
% and J Poleszczuk 2018, NCT Heidelberg
% jakob.kather?gmail.com
% 
clear variables; %clearing variables from the workspace
addpath('Subroutines'); %add necessary paths

currDname = 'timelines_KNT_CEA_20'; %name of the file with data
data = load(['./Data files/',currDname,'.mat']); % load raw data for fitting
% options for optimization
cnst.maxOpt = 1000; % default 1E3 on desktop, 1E2 on mobile

% optional: jitter
jitter.level = 0.1; % default 0.01
jitter.n = 25; % default 50 on desktop, 10 on mobile
rng('shuffle');

for jit = 1:jitter.n
data = jitterData(data,jitter); % add jitter to data
disp(['*****',newline,'starting iteration ',num2str(jit)]);

%get initial guess for the optimization
[initGuess,lb,~] = getInitGuess(data.patientCollection);
ub = []; % [] for no upper bound

%defining options for the optimizer
% WAS 1000 iterations
optimizationOpts = optimoptions('lsqnonlin','MaxIterations',cnst.maxOpt,... % optimizer options
    'MaxFunctionEvaluations',1e6,'Display','iter','Diagnostics','off',...
    'UseParallel',true, 'DiffMinChange',1e-4); 

%reformat data
patients = struct();
for i = 1:numel(data.patientCollection) %for each patient
    patients(i).timelineBIO{1} = data.patientCollection{i};
    patients(i).newID = strrep(strrep(strrep(data.patientNames{i},'LE_CRC_',''),'_BIOCHEM',''),'_','-');
end

%% STEP 1
%run the global optimizer for the first time
B = lsqnonlin(@lossfun, ones(size(initGuess)), lb, [], optimizationOpts, patients, initGuess);
save(['./Tmp/',currDname,'.mat'],'B','initGuess','patients');

%% STEP 2
%for each patient try to improve the fit while holding global (not patient specitfic) parameters
Bindividual = cell(1,numel(patients));
parfor i = 1:numel(patients) %for each patient
    %defining options for the optimizer
    optimizationOptsLoc = optimoptions('lsqnonlin','MaxIterations',cnst.maxOpt,... % optimizer options
        'MaxFunctionEvaluations',1e6,'Display','none','Diagnostics','off',...
        'UseParallel',false, 'DiffMinChange',1e-4); 

    disp(['Refining fit for patient ' patients(i).newID]);
    %check if patient has enough data points for refinement
    indx = (1:5) + (i-1)*5 + 2; %indices for patient-specific parameters 
    if size(patients(i).timelineBIO{1},1) >  4
        globalParams = B([1 2]).*initGuess([1 2]);
        x0 = B(indx); %global parameters are the first two one in B
        initGuessIndividual = initGuess(indx);
        Bindividual{i} = lsqnonlin(@lossfun, x0, lb(indx), [], optimizationOptsLoc, patients(i), initGuessIndividual, globalParams);
    else
        disp('Not enough data points for refinement')
        Bindividual{i} = B(indx); %take the same values as before
    end
end
save(['./Tmp/',currDname,'Individual.mat'],'Bindividual');

%% STEP 3
%run the global optimizer with refined initial guesses
%first we need to prepare new starting point
startingPoint = B([1 2]); %take global params
for i = 1:numel(patients) %for each patient
   startingPoint = [startingPoint; Bindividual{i}]; 
end
[Bfinal,resnorm] = lsqnonlin(@lossfun, startingPoint, lb, [], optimizationOpts, patients, initGuess);

disp('finished this iteration');
allBs(:,jit) = Bfinal;  % save the fitted parameters
allIs(:,jit) = initGuess; % save the initial guesses (should be identical)
allQs(:,jit) = resnorm; % save the goodness of fit
disp('saved data');
end

save(['./Fit results/',currDname,'_jitter.mat'],...
    'Bfinal','allBs','allQs','allIs','initGuess','patients');
