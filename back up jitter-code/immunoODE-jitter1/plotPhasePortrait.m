% JN Kather 2018, Moffitt Cancer Center, Tampa FL
% and J Poleszczuk 2018, NCT Heidelberg
% jakob.kather?gmail.com
% 
clear all, close all, clc   % clean up before run
addpath(genpath('./subroutines'));        % add subroutines
currDname = 'timelines_KNT_CEA_20'; % %'timelines_Le_CEA_30' timelines_NCT_CEA_24
load(['./Fit results/',currDname,'.mat']);                % will load previously saved data (res and params)
%initCond = ones(size(initCond)); % overwrite initial conditions
[B,params] = rescaleParamVector(Bfinal,initGuess,numel(patients));
params = params';
sq = @(varargin) varargin'; % define an auxiliary function to process cells
% first, define the parameter names
paramnames = {'\alpha (T growth)','\beta (T capacity)','\gamma (killing)','\delta (E influx)','\epsilon (E saturation)','\theta (exhaustion)','\zeta (E decay)','T0','E0'};
allIDs = sq(patients.newID); % extract patient IDs
% now, some more preferences
tmax = 180;         % how many days for the trajectory
EPlotMax = 3;
TPlotMax = 3;
Npoints = 7;        % how many arrows in the quiver plot for each dimension
numLines = 5; % default 20
lineWid = 3; % default 3
plotAlpha = 0.5;%0.01; % default 0.01
%disturb = 1+(rand(nLines,2)-0.5)*2; % how much to wiggle the trajectories
dosave = 0;         % save result as hi-res PNG image?
rng('default'); 
figure
for pat = 1:size(params,1)
   
subplot(4,ceil(size(params,1)/4),pat);

%%defining vector vield
pv = params(pat,:);
rhs = @(t,x)([pv(1).*x(1,:).*(1-pv(2).*x(1,:)) - pv(3).*(x(1,:).*x(2,:)); ...
             (pv(4).*x(1,:).*x(2,:))./(pv(5)+x(1,:)) - (pv(6).*x(1,:).* x(2,:)) - pv(7) .*x(2,:)]);

x = linspace(0,TPlotMax,Npoints);
y = linspace(0,EPlotMax,Npoints);

dx = x(2)-x(1);
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

%%calculating steady states
SS = calculateSteadyStates(pv);
stab = calculateStability(SS,pv);

%%evaluating separating curve
%yet to be done

%%Evaluating and plotting solutions
hold on
% FUNCTION RETURNING MODEL SOLUTION ON [0,tmax] FOR GIVEN INITIAL CONDITION
%solve = @(init)(ode45(rhs,[0 tmax],init));   
[XL, YL] = meshgrid(linspace(0.1,TPlotMax,numLines),linspace(0.1,EPlotMax,numLines));
CurrInitCond = [XL(:),YL(:)];%disturb.*initCond(pat,:); % get acutal initial conditions
sols = cell(1,size(CurrInitCond,1));
disp('plotting all trajectories...');
parfor i = 1:size(CurrInitCond,1) % use PARFOR
    tic
    disp(['starting iteration ',num2str(i),' of ',...
        num2str(size(CurrInitCond,1)),' at initial cond ',num2str(CurrInitCond(i,:))]);
    SSloc = SS;
    sols{i} = solveModelForPortraitWall(SSloc(stab,:), pv, CurrInitCond(i,:));%solve(CurrInitCond(i,:));
    disp(['finished iteration ',num2str(i),' of ',num2str(size(CurrInitCond,1))]);
    toc
end
disp('done with all trajectories');

% generate the colors based on the location of the stable point
% mycmap = flipud(lines(2));
% if numel(stab)==3
%     mycmap(1,:) = mycmap(2,:);
% end

stabSS = SS(stab,1);

mycmap = zeros(numel(stabSS),3); % this is stable disease

diseaseSD = ((stabSS<1.2) & (stabSS>0.7));
diseasePR = stabSS<0.7;
diseasePD = stabSS>1.2;

mycmap(diseaseSD,:) = repmat([255 166 41]/255,sum(diseaseSD),1);
mycmap(diseasePR,:) = repmat([0 0 230]/255,sum(diseasePR),1);
mycmap(diseasePD,:) = repmat([194 0 58]/255,sum(diseasePD),1);

hold on

%%plotting each solution
for i = 1:length(sols) 
    p = plot(sols{i}.y(1,:),sols{i}.y(2,:),'k-','LineWidth',lineWid);
    try
    if ~isempty(sols{i}.ie) %if it converged
        p.Color = [mycmap(sols{i}.ie,:) plotAlpha];
    else
        p.Color = [.5 .5 .5 plotAlpha];  
    end
    catch
        warning('error setting color');
        p.Color = [.5 .5 .5 plotAlpha];  
    end
end


%%plotting vector field
[X1, Y1] = meshgrid(linspace(0,TPlotMax,Npoints),linspace(0,EPlotMax,Npoints));
q = quiver(X1,Y1,U,V); %plotting vector field
q.AutoScaleFactor = 0.4; % 0.75
q.LineWidth = .8;
q.MaxHeadSize = 500 * q.MaxHeadSize;
q.Color = [.2 .2 .2];

%%plotting steady states - truncated to plotting area
ssSize = 40;
scatter(min(SS(stab,1),TPlotMax),min(SS(stab,2),EPlotMax),ssSize,'kd','filled');
%scatter(min(SS(~stab,1),TPlotMax),min(SS(~stab,2),EPlotMax),ssSize,'ko','filled');
scatter(1,1,ssSize*2,'kx','LineWidth',2)

% plot the actual solution from starting point [1,1]
disp('plotting actual solution');
solActual = solveModelForPortraitWall(SS(stab,:), pv, [1,1],1E-4); %solve(CurrInitCond(i,:));
p = plot(solActual.y(1,:),solActual.y(2,:),'k-','LineWidth',1.5);
disp('done');

% plot the solution at 60 days, pv exists in the workspace
rhs = @(t,x)([pv(1).*x(1,:).*(1-pv(2).*x(1,:)) - pv(3).*(x(1,:).*x(2,:)); ...
             (pv(4).*x(1,:).*x(2,:))./(pv(5)+x(1,:)) - (pv(6).*x(1,:).* x(2,:)) - pv(7) .*x(2,:)]);
%for tplot = (1:7)*52 %(1:4)*30      
%sol60 = ode45(rhs,[0,tplot],[1,1]);
%scatter(sol60.y(1,end),sol60.y(2,end),ssSize*0.75,'ko','filled','LineWidth',0.5);
%scatter(sol60.y(1,end),sol60.y(2,end),ssSize,'w.','LineWidth',1.5);
%end

%axis square tight
%axis([0 max(max(SS(:,1))*1.1,2) 0 max(max(SS(:,2))*1.1,2)]);
axis([0 TPlotMax 0 EPlotMax]);
axis square 
set(gcf,'Color','w');
xlabel('TUM');ylabel('EFF');
title(allIDs{pat});
drawnow

end
suptitle('probability of phase trajectories');

if dosave % save resulting figure as PNG
  set(gcf,'Position',[ 650         209        1379         882]);
print(gcf,['./output/params_phase_aug_',currDname,'.png'],'-dpng','-r400');
print(gcf,['./output/params_phase_aug_',currDname,'.svg'],'-dsvg');
end