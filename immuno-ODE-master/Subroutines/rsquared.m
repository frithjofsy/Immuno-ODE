function R2 = rsquared(yData,yFit)


    SSres = sum((yData(:)-yFit(:)).^2);
    SStot = sum((yData(:)-mean(yData(:))).^2);
    R2 = (1 - (SSres/SStot)); 
    
% yData
% yFit
% SSres
% SStot
% R2
% disp('##');

end