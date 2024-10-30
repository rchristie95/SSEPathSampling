function idx = ShootPointRND(path)
    % This function selects a random shooting point along the given trajectory path.
    % Inputs:
    %   path - The given trajectory path
    % Output:
    %   idx  - Index of the shooting point selected

    trajectorySteps = length(path);
    idx = randi([2, trajectorySteps - 1]);
end
