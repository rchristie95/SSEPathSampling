function [adjusted_weights, unbiasedhistograms] = adjust_histogram_weightsTPS(histograms, binEdges, overlaps)
    % This function adjusts histogram weights considering overlaps and calculates unbiased histograms.
    % Inputs:
    %   histograms - Matrix containing histograms for different intervals
    %   binEdges   - Array containing the edges of histogram bins
    %   overlaps   - Array specifying overlapping intervals for adjustments
    % Outputs:
    %   adjusted_weights  - Adjusted weights for overlapping intervals
    %   unbiasedhistograms - Final unbiased histogram obtained after adjustment

    % Get the number of bins and overlaps
    nBins = size(histograms, 2);
    nOverlaps = size(overlaps, 2);

    % Initialize a cell array to store the indices for each pair of bounds
    indicesCellArray = cell(1, nOverlaps);
    histograms2 = histograms;
    binCenters = binEdges(1:end-1) + 0.5 * (binEdges(2) - binEdges(1));

    % Find indices between overlap bounds for each interval
    for n = 1:nOverlaps
        lowerBound = overlaps(1, n);
        upperBound = overlaps(2, n);
        indicesCellArray{n} = findIndicesBetweenValues(binEdges, lowerBound, upperBound);
        firstNZ(n) = indicesCellArray{n}(1);
        lastNZ(n) = indicesCellArray{n}(end);
    end

    % Forward adjustment of histogram weights based on overlap means
    for n = 1:nOverlaps
        adjusted_weights(n) = mean(histograms2(n, indicesCellArray{n})) / mean(histograms2(n+1, indicesCellArray{n}));
        histograms2(n+1, :) = adjusted_weights(n) * histograms2(n+1, :);
        histograms2(n+1, 1:firstNZ(n)-1) = histograms2(n, 1:firstNZ(n)-1);
    end

    % Backward adjustment for bins beyond overlap regions
    for n = nOverlaps:-1:1
        scaleFactor = mean(histograms2(n+1, indicesCellArray{n})) / mean(histograms2(n, indicesCellArray{n}));
        histograms2(n, lastNZ(n)+1:end) = scaleFactor * histograms2(n+1, lastNZ(n)+1:end);
    end

    % Compute the unbiased histogram as the mean of all adjusted histograms
    unbiasedhistograms = mean(histograms2, 1);
    unbiasedhistograms = unbiasedhistograms / trapz(binCenters, unbiasedhistograms);

end

% Helper function to find indices between two values
function indices = findIndicesBetweenValues(array, lowerBound, upperBound)
    % Find the indices of an array between two values
    %
    % Inputs:
    %   array      - The input array
    %   lowerBound - The lower bound of the desired range
    %   upperBound - The upper bound of the desired range
    % Outputs:
    %   indices - A vector containing the indices of the array elements that fall within the specified range

    % Find the logical indices of the elements between the bounds
    logicalIndices = (array >= lowerBound) & (array < upperBound);

    % Convert logical indices to actual indices
    indices = find(logicalIndices);
end
