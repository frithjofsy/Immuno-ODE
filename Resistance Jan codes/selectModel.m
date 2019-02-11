function model = selectModel(whichModel)


    switch whichModel
        
        case 0
            
            model.solve = @(params, tmesh)((1-params.f)*exp(-params.d*tmesh)+params.f*exp(params.g*tmesh));
            model.whichFixed = {'d'};
            model.whichParamsFit = {'f','g'};
            model.nominalFixed = 0.43; %order matters
            model.nominalFitted = [0.2, 1]; %order matters
            
        case 1
            
            model.solve = @(params, tmesh)(solveModel(params,tmesh));
            model.whichParamsFit = {};
            model.whichFixed = {'K','d2'};
            model.whichParamsFit = {'f','r1','r2','d1'};
            model.nominalFixed = [10, 0]; %order matters
            model.nominalFitted = [0.2,0.92, 0.5, 0.2]; %order matters
    end
    
    function sol = solveModel(params, tmesh)
        
        sol = ode45(@odes, [0 tmesh(end)],[1-params.f; params.f],[],params);
        sol = deval(sol, tmesh);
        sol = sum(sol);
        
    end
    
    function dy = odes(~,x,params)
                dy = zeros(2,1);
                T = sum(x);
                dy(1) = params.r1*x(1)*(1 - T/params.K) - params.d1*x(1);
                dy(2) = params.r2*x(2)*(1 - T/params.K) - params.d2*x(2);
    end
            

end

