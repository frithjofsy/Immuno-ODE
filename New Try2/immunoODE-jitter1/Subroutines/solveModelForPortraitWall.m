function sol = solveModelForPortraitWall(SS, pv, init,varargin)


   %calculate minimal distance to SS now
   ID = sqrt(sum((init(:)'-SS).^2,2));
   mD = min(ID);
   defaultMinDist = 0.5; % default 0.1
   t1 = clock();
   maxCompTime = 5; %maximal computation time in seconds
   
   %iD =  (init(1)-SS(stable,1)).^2;

    rhs = @(t,x)([pv(1).*x(1,:).*(1-pv(2).*x(1,:)) - pv(3).*(x(1,:).*x(2,:)); ...
             (pv(4).*x(1,:).*x(2,:))./(pv(5)+x(1,:)) - (pv(6).*x(1,:).* x(2,:)) - pv(7) .*x(2,:)]);

    tic
    function [value,isterminal,direction] = event(~,y)
        t2 = clock(); 
        currTime = etime(t2,t1) - maxCompTime;
        value = sqrt(sum((y(:)'-SS).^2,2)) - min(defaultMinDist,mD/2);
        isterminal = ones(size(value));
        direction = 0;
        if currTime>0
            value = 0*value;
            isterminal = ones(size(value));
        end
    end

    opts = odeset('Events',@event);


    if nargin>3
        opts.RelTol = varargin{1};
        opts.AbsTol = varargin{1};
    else
        opts.AbsTol = 1E-4;
        opts.RelTol = 1E-4;
    end
    
    sol = ode45(rhs,[0 1e5],init, opts); % was [0 1e5]
    
end

