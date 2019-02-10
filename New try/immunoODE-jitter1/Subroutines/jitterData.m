function data = jitterData(data,jitter)

% create backup of original data if there is no backup so far
if ~isfield(data,'backup')
    data.backup = data.patientCollection;
end

data.patientCollection = []; % remove old data
% use the backup to jitter data
for i = 1:numel(data.backup)
    currData = data.backup{i};
    % do the jitter
    currData(:,2) = currData(:,2) ...
        .* (1 + jitter.level * randn(1,size(currData,1)))';
    data.patientCollection{i} = currData;
end

end