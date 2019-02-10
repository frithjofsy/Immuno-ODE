function SS = calculateSteadyStates( pv )

% T0 = 10^8;
% E0 = 10^6;
% 
% %s = 0;
% 
% par.sigma = 0;%s/(n*params.N0*params.N0);
% par.rho = pv(4)/(pv(3)*E0); %there is additional scaling
% par.eta = pv(5); %no additional scaling %params.g/params.N0;
% par.mu = pv(6)/(pv(3)*E0);
% par.delta = pv(7)/(pv(3)*E0);
% par.alpha = pv(1)/(pv(3)*E0);
% par.beta = pv(2);

%par.omega = params.q/(params.N0*n); %assuming E0 = N0 = D0
%par.gamma = params.n/(n*params.N0);

%defining polynomial coefficients
%C3 = par.mu*par.beta;
%C2 = -par.mu+par.beta*(par.eta*par.mu+par.delta-par.rho);
%C1 = par.sigma/par.alpha+par.rho-par.eta*par.mu-...
%    par.delta+par.beta*par.delta*par.eta;
%C0 = par.eta*(par.sigma/par.alpha-par.delta);

C2 = pv(6);
C1 = pv(6)*pv(5)-pv(4)+pv(7);
C0 = pv(7)*pv(5);

%defining polynomial to be solved
%if C3~=0
%    C = [C3 C2 C1 C0];
%else
    C = [C2 C1 C0];
%end
yst = roots(C);
yst(imag(yst)~=0) = []; %deleting complex roots
yst = real(yst);
%yst

yst(yst<0) = []; %deleting negative roots
%evaluating second coordinate of stationary state
xst = pv(1)/pv(3)*(1-pv(2)*yst); %this is steady state for effector cells, this is ok

%deleting negative roots is
yst(xst<0) = [];
xst(xst<0) = [];

SS = [yst(:) xst(:)];

SS = [0 0; SS]; %trivial steady state
SS = [1/pv(2) 0; SS];%semi trivial steady state


end

