function LangevinTPSSymmetry(h2, h4, Gamma, kT, LambdaVals, windowMat, A, n_equilib, nPathsPerWindow)
    % This function performs Langevin dynamics using two-point sampling symmetry for umbrella and unconstrained sampling.
    % Inputs:
    %   h2, h4        - Potential coefficients for quadratic and quartic terms
    %   Gamma         - Damping coefficient
    %   kT            - Thermal energy (Boltzmann constant * Temperature)
    %   LambdaVals    - Array of lambda values for umbrella sampling windows
    %   windowMat     - Matrix specifying window boundaries for umbrella sampling
    %   A             - Boundaries defining the left state
    %   n_equilib     - Number of equilibration paths per window
    %   nPathsPerWindow - Number of paths per window for sampling

    % Parameters
    B = flip(abs(A)); % Boundaries for the right state
    Xmin = -sqrt(h2 / (2 * h4)); % Minimum value of X
    V = @(x) (h4 * x.^4 - h2 * x.^2); % Potential function
    nWindows = size(windowMat, 2); % Number of windows for Umbrella Sampling
    ChoiceString = {'S', 'T', 'P', 'PT'};
    dt = 0.001; % Time step
    tf = 15; % Time duration of the path
    trajectorySteps = round(tf / dt); % Number of Langevin steps
    N_t_prime = trajectorySteps;
    tickrate=2; % Metropolis lag parameter

    Overlaps(1, :) = windowMat(1, 2:end);
    Overlaps(2, :) = windowMat(2, 1:end-1);
    momentumPertUmb = 0.01; % Perturbation for umbrella sampling
    momentumPertUnc = 0.1; % Perturbation for unconstrained sampling
    filenameWorkSpace = sprintf('SymDeepworkspaceLangevinUmbrella_T=%.2f_Gamma=%.2f.mat', kT, Gamma);

    %% Umbrella Sampling

    binEdges = windowMat(1, 1):0.05:windowMat(2, end);
    nBins = length(binEdges) - 1;
    binCenters = binEdges(1:end-1) + 0.5 * (binEdges(2) - binEdges(1));
    BRange = find((binCenters >= B(1)) & (binCenters <= B(2)));
    nRounds = 2; % Number of rounds per window
    nTotalIterations = nWindows * nRounds; % Total iterations in the parfor loop
    tempHistograms = cell(nTotalIterations, 1); % Cell array to store histograms

    parfor iter = 1:nTotalIterations
        count = 0;
        window = mod(iter - 1, nWindows) + 1;
        Round = floor((iter - 1) / nWindows) + 1;
        accepted = 0;
        stuck = 0;
        path = linspace(Xmin, LambdaVals(window), trajectorySteps);
        Psi = [path; kT * randn(1, trajectorySteps)];
        converged = false;

        while ~converged
            idx = ShootPointRND(path); % Choose shooting point
            eta = randn(1, trajectorySteps);
            Zin = [Psi(1, idx); Psi(2, idx) + momentumPertUmb * randn];
            [newPath, newPsi] = LangevinDynamicsTPS(Zin, idx, dt, tf, eta, Gamma, kT, h4, h2);
            if newPath(1) < A(2) && newPath(1) > A(1) && windowMat(1, window) < newPath(trajectorySteps) && newPath(trajectorySteps) < windowMat(2, window)
                converged = true;
                path = newPath;
                Psi = newPsi;
            end
        end

        localHistogram = zeros(1, nBins);
        while accepted < nPathsPerWindow + n_equilib
            count = count + 1;
            PaccString=sprintf('Repeat');
            stuck = stuck + 1;
            idx = ShootPointRND(path); % Choose shooting point
            eta = randn(1, trajectorySteps);
            Zin = [Psi(1, idx); Psi(2, idx) + momentumPertUmb * randn];
            [newPath, newPsi] = LangevinDynamicsTPS(Zin, idx, dt, tf, eta, Gamma, kT, h4, h2);
            [choice, Ns] = PathCheckerUmbrella(newPath, A, windowMat(:, window));
            if choice ~= 0
                [newPath, newPsi] = PathTransformerLangevin(newPath, newPsi, choice);
                [revPath, ~] = PathTransformerLangevin(path, Psi, choice);
                [~, barNs] = PathCheckerUmbrella(revPath, A, windowMat(:, window));
                Pac = PaccLangevinTransformed(newPsi, Psi, idx, kT, Gamma, dt, h4, h2, choice);
                Pac = (Ns / barNs) * Pac;
                PaccString = sprintf('Pacc_%s=%.2g', ChoiceString{choice}, Pac); % Update acceptance string
                if rand < Pac
                    stuck = 0;
                    accepted = accepted + 1;
                    path = newPath;
                    Psi = newPsi;
                end
            end
            if accepted > n_equilib && mod(count, tickrate) == 0
                localHistogram = localHistogram + histcounts(path(end), binEdges);
            end
            fprintf('Langevin window %d, round %d, Progress: %5.2f, Stuck: %5.d, %s\n', window, Round, (accepted - n_equilib) / nPathsPerWindow, stuck, PaccString);
        end
        fprintf('Langevin window %d, round %d finished with accept ratio: %.2f\n', window, Round, accepted / count);
        tempHistograms{iter} = localHistogram;
    end

    % Combine histograms from cell array after parfor loop
    histograms = zeros(nWindows, length(binCenters));
    for iter = 1:nTotalIterations
        window = mod(iter - 1, nWindows) + 1;
        histograms(window, :) = histograms(window, :) + tempHistograms{iter};
    end

    normalized_histograms = zeros(nWindows, nBins);
    for window = 1:nWindows
        normalized_histograms(window, :) = histograms(window, :) / trapz(binCenters, histograms(window, :));
    end
    save(filenameWorkSpace);

    [~, unbiased_histogram_Langevin] = adjust_histogram_weightsTPS(normalized_histograms, binEdges, Overlaps);
    Ct_Umbrella_Langevin = trapz(binCenters(BRange), unbiased_histogram_Langevin(BRange));
    unbiased_free_energy_Langevin = -kT * log(unbiased_histogram_Langevin);
    unbiased_free_energy_Langevin = real(unbiased_free_energy_Langevin - min(unbiased_free_energy_Langevin));
    save(filenameWorkSpace);
    disp("Langevin umbrella sampling time below");
    toc

    %% Unconstrained Sampling

    nPathsPerWindow = 3*nPathsPerWindow;
    filenameWorkSpace = sprintf('SymDeepworkspaceLangevinUnconstrained_T=%.2f_Gamma=%.2f.mat', kT, Gamma);
    tf = 40; % Time duration of the path
    trajectorySteps = round(tf / dt); % Number of Langevin steps
    nCopy = 16;
    acceptRatio = zeros(1, nCopy);
    numGroups = 200;
    totalPaths = zeros(1, nCopy);

    tic
    h_AB_win = zeros(nCopy, numGroups);
    parfor copy = 1:nCopy
        count = 0;
        accepted = 0;
        path = zeros(1, trajectorySteps);
        Psi = zeros(2, trajectorySteps);
        stuck = 0;
        while path(1) > A(2) || ~any(path > B(1))
            eta = randn(1, trajectorySteps);
            [path, Psi] = LangevinDynamicsTPS([0; sqrt(kT) * randn], round(trajectorySteps / 2), dt, tf, eta, Gamma, kT, h4, h2);
        end
        while accepted < nPathsPerWindow + n_equilib
            count = count + 1;
            PaccString=sprintf('Repeat');
            stuck = stuck + 1;
            idx = ShootPointRND(path); % Choose shooting point
            eta = randn(1, trajectorySteps);
            Zin = [Psi(1, idx); Psi(2, idx) + momentumPertUnc * randn];
            [newPath, newPsi] = LangevinDynamicsTPS(Zin, idx, dt, tf, eta, Gamma, kT, h4, h2);
            [choice, Ns] = PathCheckerUnconstrained(newPath, A, B);
            if choice ~= 0
                [newPath, newPsi] = PathTransformerLangevin(newPath, newPsi, choice);
                [revPath, ~] = PathTransformerLangevin(path, Psi, choice);
                [~, barNs] = PathCheckerUnconstrained(revPath, A, B);
                Pac = PaccLangevinTransformed(newPsi, Psi, idx, kT, Gamma, dt, h4, h2, choice);
                Pac = (Ns / barNs) * Pac;
                PaccString = sprintf('Pacc_%s=%.2g', ChoiceString{choice}, Pac);
                if rand < Pac
                    accepted = accepted + 1;
                    path = newPath;
                    Psi = newPsi;
                    stuck = 0;
                end
            end
            if accepted > n_equilib && mod(count, tickrate) == 0
                totalPaths(copy) = totalPaths(copy) + 1;
                h_AB_win(copy, :) = h_AB_win(copy, :) + mean(reshape(indicatorB(path, B(2), B(1)), [], numGroups), 1);
            end
            fprintf('Langevin copy %d, Progress: %5.2f, Stuck: %5.d, %s\n', copy , (accepted - n_equilib) / nPathsPerWindow, stuck, PaccString);
        end
        h_AB_win(copy, :) = h_AB_win(copy, :) / totalPaths(copy);
        acceptRatio(copy) = accepted / totalPaths(copy);
        fprintf('Langevin copy %d, finished with accept ratio: %.2f\n', copy, acceptRatio(copy));
    end
    fprintf('Langevin unconstrained mean accept ratio: %.2g\n', mean(acceptRatio));

    disp("Langevin unconstrained sampling time below");
    toc
    save(filenameWorkSpace);

    %% Processing

    h_Star_AB_Langevin = mean(h_AB_win, 1);
    N_t_smooth = N_t_prime * (numGroups / trajectorySteps);
    smoothSpan = linspace(0, tf, numGroups);
    dh_Star_AB_Langevin = central_diff_8th_order(smoothSpan, h_Star_AB_Langevin);
    Correl_Langevin = h_Star_AB_Langevin / (h_Star_AB_Langevin(N_t_smooth)) * Ct_Umbrella_Langevin;
    rate_Langevin = dh_Star_AB_Langevin / (h_Star_AB_Langevin(N_t_smooth)) * Ct_Umbrella_Langevin;
    save(filenameWorkSpace);
end
