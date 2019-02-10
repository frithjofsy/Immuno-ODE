function [B,params] = rescaleParamVector(Bin,initGuess,numPat)

B = Bin.*initGuess; % rescale parameters. remember: the optimizer will
% not work with the actual vector of parameters but with a vector of
% scaling factors acting on the initial guesses. The point of this is
% to bring all parameters to the same range.
% store params
params = reshape(B((3):end),[],numPat);
params = [repmat(B(1),1,size(params,2));  params; repmat(B(2),1,size(params,2))];

end