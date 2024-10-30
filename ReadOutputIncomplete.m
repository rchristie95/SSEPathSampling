function [RateLangevinLinear,RateSSELinear]=ReadOutputIncomplete(filename,TempRange,GammaRange,nAvg)
nSweepT=length(TempRange);
nSweepG=length(GammaRange);

for n=1:nSweepT
    TempRangeRound(n)=str2double(sprintf('%.4g', TempRange(n)));
end

for n=1:nSweepG
    GammaRangeRound(n)=str2double(sprintf('%.4g', GammaRange(n)));
end

fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open file %s', filename);
end

% Read the file line by line
lines = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
lines = lines{1};

% Extract the lines containing the array data
dataLines = lines(~cellfun(@isempty, regexp(lines, 'T_rate')));

% Parse the data
SSECount=0;
LangevinCount=0;

SSECount = 0;
LangevinCount = 0;
for i = 1:length(dataLines)
%     tokens = regexp(dataLines{i}, 'T_rate_(\w+):([\d.]+),\s*kT:([\d.]+),\s*Gamma:([\d.]+)', 'tokens');
    tokens = regexp(dataLines{i}, 'T_rate_(\w+):([\de.-]+),\s*kT:([\de.-]+),\s*Gamma:([\de.-]+)', 'tokens');
    if isempty(tokens)  % add this line to handle the case where no match is found
        continue;
    end
    tokens = tokens{1};
    if strcmp(tokens{1}, 'SSE')
        % Store the data in the array
        % Column 1: T_rate, Column 2: kT, Column 3: Gamma
        SSECount = SSECount + 1;
        SSEArray(SSECount, 1) = str2double(tokens{2});
        SSEArray(SSECount, 2) = str2double(tokens{3});
        SSEArray(SSECount, 3) = str2double(tokens{4});
    elseif strcmp(tokens{1}, 'Langevin')
        % Store the data in the array
        % Column 1: T_rate, Column 2: kT, Column 3: Gamma
        LangevinCount = LangevinCount + 1;
        LangevinArray(LangevinCount, 1) = str2double(tokens{2});
        LangevinArray(LangevinCount, 2) = str2double(tokens{3});
        LangevinArray(LangevinCount, 3) = str2double(tokens{4});
    end
end

SSEArraywoReplacement=SSEArray;
LangevinArraywoReplacement=LangevinArray;

totalIterations=nSweepT*nSweepG*nAvg;
RateSSELinear = nan(totalIterations, 1);
RateLangevinLinear = nan(totalIterations, 1);

for idx=1:totalIterations
    trueIdx = ceil(idx / nAvg);
    m = mod(trueIdx-1, nSweepT) + 1;
    k = ceil(trueIdx / nSweepT);
    Gamma = GammaRangeRound(k);
    kT = TempRangeRound(m);
    logicalIndex = (LangevinArraywoReplacement(:, 2) == kT) & (LangevinArraywoReplacement(:, 3) == Gamma);
    firstMatchIndex = find(logicalIndex, 1);
    if ~isempty(firstMatchIndex)
        RateLangevinLinear(idx)=LangevinArraywoReplacement(firstMatchIndex,1);
        LangevinArraywoReplacement(firstMatchIndex,:)=[];
    end
    logicalIndex = (SSEArraywoReplacement(:, 2) == kT) & (SSEArraywoReplacement(:, 3) == Gamma);
    firstMatchIndex = find(logicalIndex, 1);
    if ~isempty(firstMatchIndex)
        RateSSELinear(idx)=SSEArraywoReplacement(firstMatchIndex,1);
        SSEArraywoReplacement(firstMatchIndex,:)=[];
    end
end

end