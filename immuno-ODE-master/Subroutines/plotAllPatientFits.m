function [initCond,BCresponses] =...
    plotAllPatientFits(res,solutions,cnst)
    
    figure(102)
    
    T0 = 1; % we are working with relative changes of cell numbers, so the 
    E0 = ones(numel(res),1);
    
    for pat = 1:numel(res) % for each patient
       
       t1 = res(pat).timelineBIO{1}(:,1); % time points for biomarker measurement
       y1 = res(pat).timelineBIO{1}(:,2)/100+1; % convert percent to relative
       subplot(cnst.numRow,ceil(numel(res)/cnst.numRow),pat) % start plot for this patient
       hold on
       % plot data and model fit
       if ~isempty(solutions{pat})
       xi = solutions{pat}.x; % retrieve time points 
       Ti = T0*solutions{pat}.y(1,:); % retrieve solutions for tumor
       Ei = E0(pat)*solutions{pat}.y(2,:); % retrieve solutions for effector
       plot(xi, Ei,'LineWidth',1,'Color','b'); % plot effector timeline
       scatter(t1,T0*y1,25,'kd','filled'); % plot the actual tumordata 
       plot(xi, Ti,'LineWidth',1,'Color','r'); % plot tumor timeline on top
       
       plot(xi, Ei+Ti,'LineWidth',2.5,'Color','k'); % plot tumor timeline on top
       
       
         % add decorations
       hold off
       set(gca,'FontSize',10);
       axis square tight
       xlabel('time (d)'); ylabel('\Delta cells');
       
        initCond(pat,1:2) = [T0,E0(pat)]; % save initial conditions 
        BCresponses(pat) = y1(end)/y1(1); % save biochemical responses
        % calculate R2 (coefficient of determination)
        ymodelERR = deval(solutions{pat}, t1);% evaluate at time points t1
        ymodelERR = ymodelERR(1,:);
        rsq(pat) = rsquared(y1,ymodelERR);
        descr1 = R2string(rsq(pat),'T'); % TUM goodness of fit
        
        % decorations
        title([res(pat).newID ,descr1]); %,newline
        xlim([0 max(t1)+1])
        ylim([0 max(T0*y1)*1.5]);
        
       else
           warning(['there was no solution for patient ',num2str(pat)]);
       end
               
        % draw line at N days if desired
        if isfield(cnst,'plotVerticalLine')
            hold on
            plot([cnst.plotVerticalLine,cnst.plotVerticalLine],...
                [0,max(T0*y1)*1.5],'k-','LineWidth',1.2);
            hold off
        end
       
        drawnow
    end
    
    
%     figure(102)
%     try %workaround for matlab graphics error if one panel is empty
%     subplot(cnst.numRow,ceil(numel(res)/cnst.numRow),i+1)
%     axis off
%     end
     set(gcf,'Color','w');
    suptitle(['mean R² ',num2str(round(mean(rsq),2)),', median R² ',num2str(round(median(rsq),2))]);
    
end