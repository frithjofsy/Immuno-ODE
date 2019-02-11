function plotResults(results)

for i = 1:numel(results)
    
   figure(i)
   clf
   hold on
   plot(results{i}.timepoints, results{i}.data, 'Marker','o','LineStyle','none');
   plot(results{i}.tmesh, results{i}.solution);
   
   hold off
   
end


end

