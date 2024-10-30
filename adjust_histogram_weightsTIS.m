function [adjusted_weights, unbiasedhistogram] = adjust_histogram_weightsTIS(histograms, XVals, binCenters, successfrac)
    % This function adjusts the histogram weights using optimization and calculates an unbiased histogram.
    % Inputs:
    %   histograms   - Matrix containing histograms for different intervals
    %   XVals        - Array containing the x-values representing interval boundaries
    %   binCenters   - Array containing the center values of histogram bins
    %   successfrac  - Success fraction used for cumulative adjustments
    % Outputs:
    %   adjusted_weights  - The optimized adjusted weights
    %   unbiasedhistogram - The computed unbiased histogram

    % Number of bins and histograms
    nBins = size(histograms, 2);
    NGrams = size(histograms, 1);

    % Compute cumulative success fractions
    successfrac = cumprod(successfrac, 2);
    successfrac = mean(successfrac, 1);
    successfrac = [1 successfrac];

    % Initialize histograms for adjustments
    histograms2 = zeros(NGrams, nBins);
    histograms2(1, :) = histograms(1, :);

    % Initialize a cell array to store the indices for each pair of bounds
    indicesCellArray = cell(1, NGrams - 1);

    % Iterate over each pair of histograms
    for n = 1:NGrams-1
        % Find indices where the bin centers fall within the interval bounds
        Indices = (binCenters >= XVals(n+1)) & (binCenters <= XVals(n+2));
        firstNZ = find(Indices, 1);

        % Update histograms with cumulative success fraction
        histograms2(n+1, :) = successfrac(n+1) * histograms(n+1, :);
        histograms2(n+1, 1:firstNZ-1) = histograms2(n, 1:firstNZ-1);

        % Store the indices in the cell array
        indicesCellArray{n} = Indices;
    end

    %% Global optimization for weight adjustment

    % Initial weights
    initial_weights = ones(NGrams, 1) / NGrams;

    % Create problem structure for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(weights) objective_function(weights, histograms2, indicesCellArray, NGrams), ...
        'x0', initial_weights, ...
        'Aeq', ones(1, NGrams), ...
        'beq', 1, ...
        'options', optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'MaxIterations', 1000000, 'MaxFunctionEvaluations', 10000000, ...
        'OptimalityTolerance', 1e-50, 'StepTolerance', 1e-30, 'Display', 'off'));

    % Create a GlobalSearch object
    gs = GlobalSearch('Display', 'off');

    % Run the optimization
    [adjusted_weights, ~, ~, ~, ~] = run(gs, problem);

    % Compute unbiased histogram
    unbiasedhistogram = sum(repmat(adjusted_weights, 1, nBins) .* histograms2, 1);

    % Adjust the weights using cumulative success fractions
    adjusted_weights = adjusted_weights .* successfrac(1:end-1).';
    adjusted_weights = adjusted_weights / sum(adjusted_weights);
end

% Objective function for optimization
function error_out = objective_function(weights, histograms, indicesCellArray, NGrams)
    % This function computes the error for given weights based on the difference between successive histograms.
    % Inputs:
    %   weights       - The weights for optimization
    %   histograms    - Adjusted histograms for each interval
    %   indicesCellArray - Cell array containing indices for each interval
    %   NGrams        - Number of histograms
    % Output:
    %   error_out     - The computed error value to be minimized

    error_out = 0;
    for n = 1:NGrams - 1
        % Calculate the difference between successive histograms weighted by current weights
        error_out = error_out + sum((histograms(n, :) * weights(n) - histograms(n+1, :) * weights(n+1)).^2);
    end
end
