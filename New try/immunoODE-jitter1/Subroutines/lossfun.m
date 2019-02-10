function [err, solOut] = lossfun( x, res, initGuess, globalParams)

if nargin < 4
   globalParams = []; 
end
numSpecificParams = 1; %number of patient specific params

x = initGuess.*x;

if nargout > 1
    solOut = cell(1,numel(res));
end

err = [];
for i = 1:numel(res) %iterate all patients
    % ---- COLLECT and PREPARE data
    
    % BIO = serum biomarker level of tumor cells. Assuming BIO ~ T
    t1 = res(i).timelineBIO{1}(:,1); % time points for BIO measurement
    y1 = res(i).timelineBIO{1}(:,2)/100+1; % convert percent to relative
    
    T = max(res(i).timelineBIO{1}(:,1)); %get maximal time
    
    %second we solve the model for the given set of parameters
    if ~isempty(globalParams)
        xLoc = [globalParams; x((1:numSpecificParams)+(i-1)*numSpecificParams)];
    else
        xLoc = [x([1 2]); x((1:numSpecificParams)+2+(i-1)*numSpecificParams)];
    end
    
    sol = solveModel(xLoc, T);
    if nargout > 1
        solOut{i} = sol;
    end
    
    if isempty(sol) %solution is corrupted
        S = 20*ones(2,length(t1));
    else
        S = deval(sol,t1);
    end
    
    % err = [err; (S(1,:)'-y1)];
   % lambda = 0.0001; % median loss for 100 days training in validation cohort, 500 iterations
    % no beta removed               removed beta
    % -----------                   --------------
    % 0->0.6                        na
    % 0.001->0.68                   0.58 (but good mean)    
    % 0.01->0.64                    0.78 (but bad mean)
    % 0.05->0.79
    % 0.1->0.8                      ###
    % 0.15->0.84 (but bad mean)
    % 0.2->0.82 (but bad mean)
    % 0.25->
    % 0.5->0.25 (but good mean)
    % 0.75->-0.05 (but bad mean)
    %
    % 0 ->-230|0.6
    % 0.1->bad|0.78
    % 0.001->bad|0.78
    % 0.0001->
    % next step: disturb input dat,a fit repeatedly, average parameters
   % patientSpecParams = x((1:numSpecificParams)+(i-1)*numSpecificParams);
    currloss = (S(1,:)'-y1);% +  sum((lambda *patientSpecParams));
    err = [err; currloss];
end

end

