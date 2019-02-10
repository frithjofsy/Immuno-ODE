function [t1,yd,xi,Ti,Ei,patnames,R2] = collectJitterFits(res,solutions)

    T0 = 1; % we are working with relative changes of cell numbers, so the 
    E0 = ones(numel(res),1);
    
    for pat = 1:numel(res) % for each patient
       
       t1{pat} = res(pat).timelineBIO{1}(:,1); % time points for biomarker measurement
       y1{pat} = res(pat).timelineBIO{1}(:,2)/100+1; % convert percent to relative
       yd{pat} = T0*y1{pat};
       % plot data and model fit
       if ~isempty(solutions{pat})
       xi{pat} = solutions{pat}.x; % retrieve time points 
       Ti{pat} = T0*solutions{pat}.y(1,:); % retrieve solutions for tumor
       Ei{pat} = E0(pat)*solutions{pat}.y(2,:); % retrieve solutions for effector

       % calculate R2 (coefficient of determination)
       ymodelERR = deval(solutions{pat}, t1{pat});% evaluate at time points t1
       ymodelERR = ymodelERR(1,:);
       R2(pat) = rsquared(y1{pat},ymodelERR);
      
       % patient name
       patnames{pat} = res(pat).newID;
       else
           xi{pat} = [];
           Ti{pat} = [];
           Ei{pat} = [];
           warning(['there was no solution for patient ',num2str(pat)]);
    end

end