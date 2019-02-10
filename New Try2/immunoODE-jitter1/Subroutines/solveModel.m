function sol = solveModel( par, Tmax )

%adding CPU wall time
maxCompTime = 5; %maximal computation time in seconds
t1 = clock();

opts = odeset('Events',@CPUwall);
sol = ode45(@myODEs, [0 Tmax], ones(2,1), opts); % use ode45 or ode23s (stiff)

if ~isempty(sol.ie) %hit CPU wall time
    sol = [];
end

    function [value, isterminal,direction] = CPUwall(~,~)
       t2 = clock(); 
       value = etime(t2,t1) - maxCompTime;
       isterminal = 1;
       direction = 1;
    end



    function dxdt = myODEs(~, x) % define the coupled ODEs
        dxdt = zeros(size(x)); % preallocate function values

        alpha = par(1); % common for each patient, previously 0.0274; 
        delta = par(2); % d common for each patient
        
        t = [1:100];
        T = x(1);
        E = x(2);
%        for i:numel(t)
        % model definition
%         dxdt(1) = (alpha*T.*(1-par(3).*T)) - par(4).*(T.*E); 
%         dxdt(2) = (par(5).*T.*E)./(par(6)+T) - (par(7).*T.* E) - delta*E;
%         dxdt(1) = (alpha*T.*(1-par(3).*T)) - par(4).*(T.*E); 
%         dxdt(2) = (par(1).*T.*E)./(par(2)+T) - (par(3).*T.* E) - delta*E;
         dxdt(1) = (1-par(3)).*par(1).*exp(par(1).*T);  %-(1-f)e^(bt)
         dxdt(2) = par(3).*par(2).*exp(par(2).*T);        % f*g*(e^gt)
%          dxdt(1) = (1-par(3)).*par(1).*exp(par(1).*t(i);  %-(1-f)e^(bt)
%          dxdt(2) = par(3).*par(2).*exp(par(2).*t(i));        % f*g*(e^gt)
          
    end
end

