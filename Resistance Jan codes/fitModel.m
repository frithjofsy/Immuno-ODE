function results = fitModel(patients, model)


    NoPatients = numel(patients);
    results = cell(1,NoPatients);
    
    for i = 1:NoPatients
       %do the fitting 
       timepoints = patients{i}.time;
       values = patients{i}.data;
       
       opts = optimset('Display','iter');
       xF = lsqnonlin(@F, model.nominalFitted, zeros(size(model.nominalFitted)), [], opts);
     
       tmesh = linspace(timepoints(1),timepoints(end),100);
       [res.err, res.solution] = F(xF);
       res.timepoints = timepoints;
       res.data = values;
       res.tmesh = tmesh;
       results{i} = res;
       
    end
    
    

    function [err, solution] = F(x)
       
        params = struct();
        
        for k = 1:length(model.whichFixed)
            params.(model.whichFixed{k}) = model.nominalFixed(k);
        end
        
        for k = 1:length(model.whichParamsFit)
            params.(model.whichParamsFit{k}) = x(k);
        end
        
        sol = model.solve(params, timepoints); 
        sol = sol(:);
        
        err = values - sol;
        
        if nargout > 1
            solution = model.solve(params, tmesh); 
        end
        
    end
    
    
end



