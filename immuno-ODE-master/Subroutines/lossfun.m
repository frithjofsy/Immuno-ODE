function [err, solOut] = lossfun( x, res, initGuess, globalParams)

if nargin < 4
   globalParams = []; 
end
% numSpecificParams = 5; %number of patient specific params
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
        xLoc = [x([1 2]); x((1:numSpecificParams)+1+(i-1)*numSpecificParams)];
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
    
  
    currloss = (S(1,:)'+S(2,:)')-y1;% +  sum((lambda *patientSpecParams));
    err = [err; currloss];
end

end

