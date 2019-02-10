function plotJitterFits(t1,yd,xi,Ti,Ei,patnames,R2,fitQ,cnst)

T0 = 1;

% t1 = time points for actual data, same for all jitter experiments
% yd = y values for actual data, same for all jitter experiments
% patnames = patient names, same for all jitter experiments
numJit = numel(xi); % number of jitter runs
numPat = numel(xi{1}); % number of patients
figure(102),clf
    for pat = 1:numPat % for each patient
       
       subplot(cnst.numRow,ceil(numPat/cnst.numRow),pat) % start plot for this patient      
       hold on
        if isfield(cnst,'plotVerticalLine')
            if strcmpi(cnst.plotVerticalLine,'half')
             plot([max(t1{1}{pat})/2,max(t1{1}{pat})/2],...
                [0,max(T0*yd{1}{pat})*1.5],'k-','LineWidth',1.2);  
            elseif strcmpi(cnst.plotVerticalLine,'HalfPoints')
             indx = ceil(numel(t1{1}{pat})/2);
             xLine = mean([t1{1}{pat}(indx),t1{1}{pat}(indx+1)]);
             plot([xLine,xLine],...
                [0,max(T0*yd{1}{pat})*1.5],'k-','LineWidth',1.2); 
            elseif strcmpi(cnst.plotVerticalLine,'crop33perc')
             indx = ceil(numel(t1{1}{pat})*(2/3));
             xLine = mean([t1{1}{pat}(indx),t1{1}{pat}(indx+1)]);
             plot([xLine,xLine],...
                [0,max(T0*yd{1}{pat})*1.5],'k-','LineWidth',1.2);       
            elseif cnst.plotVerticalLine < 0
            plot([max(t1{1}{pat})+cnst.plotVerticalLine,max(t1{1}{pat})+cnst.plotVerticalLine],...
                [0,max(T0*yd{1}{pat})*1.5],'k-','LineWidth',1.2);
            else
            plot([cnst.plotVerticalLine,cnst.plotVerticalLine],...
                [0,max(T0*yd{1}{pat})*1.5],'k-','LineWidth',1.2);
            end
        end
 
        clear allEs allTs
        
        [~,targetFit] = min(fitQ); % use the time line from the best fit
        targetX = linspace(0,max(xi{targetFit}{pat}),1000); % plot at 1000 time points
        
        for fits = 1:numJit
            try
            allEs(fits,:) = interp1(xi{fits}{pat},Ei{fits}{pat},targetX);
            allTs(fits,:) = interp1(xi{fits}{pat},Ti{fits}{pat},targetX);
            catch
            allEs(fits,:) = repmat(NaN,numel(targetX),1);
            allTs(fits,:) = repmat(NaN,numel(targetX),1);
            end
        end
      


       medianE = nanmedian(allEs,1);
       medianT = nanmedian(allTs,1);
      % stdE = nsd*std(allEs,1,'omitnan');
      % stdT = nsd*std(allTs,1,'omitnan');
      loCI_E = quantile(allEs,cnst.loQuantile);
      loCI_T = quantile(allTs,cnst.loQuantile);
      hiCI_E = quantile(allEs,1-cnst.loQuantile);
      hiCI_T = quantile(allTs,1-cnst.loQuantile);
       ysE = [loCI_E,fliplr(hiCI_E)];
       ysT = [loCI_T,fliplr(hiCI_T)];
       % crop
       ysE(ysE<0) = 0;
       ysT(ysT<0) = 0;
       ysE(ysE>max(T0*yd{1}{pat})*1.5) = max(T0*yd{1}{pat})*1.5;
       ysT(ysT>max(T0*yd{1}{pat})*1.5) = max(T0*yd{1}{pat})*1.5;
         patch([targetX,fliplr(targetX)],ysE,[0 0 1],'EdgeColor','none','FaceAlpha',0.15); % plot the confidence
         patch([targetX,fliplr(targetX)],ysT,[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

       plot(targetX,medianE,'LineWidth',2,'Color',[0 0 1 1]); % plot the median line for immune
       plot(targetX,medianT,'LineWidth',2,'Color',[1 0 0 1]); % plot the median line for tumor
           
%            % plot the fits on top
%       for fits = 1:numJit
%           alpha = 0.25;
%          plot(xi{fits}{pat}, Ei{fits}{pat},'LineWidth',0.5,'Color',[0 0 1 alpha]); % plot effector timeline
%          plot(xi{fits}{pat}, Ti{fits}{pat},'LineWidth',0.5,'Color',[1 0 0 alpha]); % plot tumor timeline on top 
%       end
    
       scatter(t1{1}{pat},T0*yd{1}{pat},25,'kd','filled'); % plot the actual tumor data 

        title(char(patnames{1}{pat})); %,newline
        xlim([0 max(t1{1}{pat})+1])
        ylim([0 max(T0*yd{1}{pat})*1.5]);
        
        hold off
        drawnow
        
    end
    
    set(gcf,'Color','w')
end