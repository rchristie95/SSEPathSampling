function SSETPSSymmetry(h2, h4, Gamma, kT, LambdaVals, windowMat, A, n_equilib, nPathsPerWindow)
    % This function performs Symmetric Sampling with Stochastic Schrödinger Equation (SSE) and TPS.
    % Inputs:
    %   h2, h4           - Potential coefficients
    %   Gamma            - Damping coefficient
    %   kT               - Thermal energy (Boltzmann constant * Temperature)
    %   LambdaVals       - Values for Lambda (Umbrella sampling parameters)
    %   windowMat        - Matrix defining window boundaries for umbrella sampling
    %   A                - Boundaries defining state A
    %   n_equilib        - Number of equilibrium paths per window
    %   nPathsPerWindow  - Number of paths per window
    %
    % Outputs:
    %   No explicit outputs, data is saved to file

    % Set up initial parameters and basis
    B = flip(abs(A)); % Define boundaries for state B, which is the flipped version of state A
    Xmin = -sqrt(h2 / (2 * h4)); % Minimum value for X, based on potential coefficients
    bSize = 40; % Basis size for the Hilbert space

    % Generate annihilation operator
    a = zeros(bSize); % Initialize the annihilation operator matrix
    for n = 1:bSize-1
        a(n, n+1) = sqrt(n); % Populate the annihilation operator with appropriate values
    end
    ChoiceString={'S','T','P','PT'}; % Symmetry options

    % Define constants and operators
    hbar = 1; % Reduced Planck's constant
    Xhat = sqrt(hbar / 2) * (a' + a); % Position operator
    Phat = 1i * sqrt(hbar / 2) * (a' - a); % Momentum operator
    Lhat = sqrt(4 * Gamma * kT / hbar) * Xhat + 1i * sqrt(Gamma * hbar / (4 * kT)) * Phat; % Lindblad operator for the system
    HhatGamma = 0.5 * Phat^2 + h4 * Xhat^4 - h2 * Xhat^2 + 0.5 * Gamma * (Xhat * Phat + Phat * Xhat); % Hamiltonian with damping term
    dt = 0.001; % Time step for integration
    tf = 15; % Time duration of the path
    trajectorySteps = round(tf / dt); % Number of SSE steps per path
    TSpan = linspace(1, tf, trajectorySteps); % Time span for plotting
    N_t_prime = trajectorySteps;
    Overlaps(1, :) = windowMat(1, 2:end); % Define overlaps for umbrella sampling
    Overlaps(2, :) = windowMat(2, 1:end-1);
    momentumPertUmb = 0.01; % Momentum perturbation for umbrella sampling
    momentumPertUnc = 0.1; % Momentum perturbation for unconstrained sampling
    tickrate=2; % Metropolis lag parameter
    ThermalState = LindbladNull(HhatGamma, Lhat, hbar); % Compute the thermal state using Lindblad dynamics
    PsiGround = zeros(bSize, 1); % Initial state (ground state)
    PsiGround(1) = 1; % Set the first element to 1 to represent the ground state

    % File name for storing workspace
    filenameWorkSpace = sprintf('SymDeepworkspaceSSEUmbrella_T=%.2f_Gamma=%.2f.mat', kT, Gamma);

    %% Umbrella Sampling
    % Define histogram bins for umbrella sampling
    binEdges = windowMat(1, 1):0.05:windowMat(2, end);
    nBins = length(binEdges) - 1;
    binCenters = binEdges(1:end-1) + 0.5 * (binEdges(2) - binEdges(1)); % Calculate bin centers
    BRange = find((binCenters >= B(1)) & (binCenters <= B(2))); % Define the range for state B
    nRounds = 2; % Number of rounds per window for sampling
    nWindows = size(windowMat, 2); % Number of windows for umbrella sampling
    nTotalIterations = nWindows * nRounds; % Total iterations in the parfor loop
    tempHistograms = cell(nTotalIterations, 1); % Using cell array to store histograms for each iteration

    parfor iter = 1:nTotalIterations
        % Initialize variables for each iteration
        count = 0;
        window = mod(iter - 1, nWindows) + 1; % Determine which window to sample
        Round = floor((iter - 1) / nWindows) + 1; % Determine the current round
        accepted = 0;
        stuck = 0;
        % Generate an initial path with added noise
        path = linspace(Xmin, LambdaVals(window), trajectorySteps) + sqrt(kT) * 1i * randn(1, trajectorySteps);
        Psi = DisplacementOpG(path / sqrt(2), a, PsiGround); % Displace the ground state along the generated path
        path = real(path); % Take only the real part of the path
        converged = false;
        % Converge to a valid path
        while ~converged
            idx = ShootPointRND(path); % Choose a random shooting point along the path
            eta = randn(1, trajectorySteps-1); % Generate random noise for SSE dynamics
            PsiPert = DisplacementOp(sqrt(kT) * 1i * randn / sqrt(2), a) * Psi(:, idx); % Apply a perturbation to the current state
            [newPath, newPsi] = SSEDynamicsTPS(PsiPert, dt, tf, eta, HhatGamma, Lhat, Xhat, hbar, idx); % Run SSE dynamics to generate a new path
            % Check if the new path satisfies the boundary conditions for state A and the current window
            if newPath(1) < A(2) && newPath(1) > A(1) && windowMat(1, window) < newPath(trajectorySteps) && newPath(trajectorySteps) < windowMat(2, window)
                converged = true; % If the conditions are met, mark convergence as true
                path = newPath;
                Psi = newPsi;
            end
        end
        localHistogram = zeros(1, nBins); % Initialize local histogram to store path data
        % Perform sampling until the required number of paths is accepted
        while accepted < nPathsPerWindow + n_equilib
            count = count + 1;
            PaccString = sprintf('Repeat'); % Default acceptance string
            stuck = stuck + 1;
            idx = ShootPointRND(path); % Choose a random shooting point along the path
            eta = randn(1, trajectorySteps-1); % Generate random noise for SSE dynamics
            PsiPert = DisplacementOp(momentumPertUmb * 1i * randn / sqrt(2), a) * Psi(:, idx); % Apply a perturbation to the current state
            [newPath, newPsi] = SSEDynamicsTPS(PsiPert, dt, tf, eta, HhatGamma, Lhat, Xhat, hbar, idx); % Run SSE dynamics to generate a new path
            [choice, Ns] = PathCheckerUmbrella(newPath, A, windowMat(:, window)); % Check if the path satisfies boundary conditions for umbrella sampling
            if choice ~= 0
                [newPath, newPsi] = PathTransformerSSE(newPath, newPsi, choice); % Transform the path based on the chosen symmetry
                [revPath, ~] = PathTransformerSSE(path, Psi, choice); % Reverse the original path for comparison
                [~, NsBar] = PathCheckerUmbrella(revPath, A, windowMat(:, window)); % Check the reverse path
                Pac = PaccSSETransformed(newPsi, Psi, idx, ThermalState, HhatGamma, Lhat, hbar, dt, choice); % Calculate the acceptance probability
                Pac = (Ns / NsBar) * Pac; % Adjust the acceptance probability based on symmetry factors
                PaccString = sprintf('Pacc_%s=%.2g', ChoiceString{choice}, Pac); % Update acceptance string
                if rand < Pac % Accept the new path with probability Pac
                    accepted = accepted + 1;
                    path = newPath;
                    Psi = newPsi;
                    stuck = 0;
                end
            end
            % Update histogram if the path is accepted after equilibrium
            if accepted > n_equilib && mod(count, tickrate) == 0
                localHistogram = localHistogram + histcounts(path(end), binEdges);
            end
            fprintf('SSE window %d, round %d, Progress: %5.2f, Stuck: %5.d, %s\n', window, Round, (accepted - n_equilib) / nPathsPerWindow, stuck, PaccString);
        end
        fprintf('SSE window %d, round %d finished with accept ratio: %.2g\n', window, Round, accepted / count);
        tempHistograms{iter} = localHistogram; % Store the local histogram in the cell array
    end
    save(filenameWorkSpace);

    % Combine histograms from all iterations after parfor loop
    histograms = zeros(nWindows, nBins);
    for iter = 1:nTotalIterations
        window = mod(iter - 1, nWindows) + 1;
        histograms(window, :) = histograms(window, :) + tempHistograms{iter};
    end

    % Normalizing histograms for each window
    for window = 1:nWindows
        normalized_histograms(window, :) = histograms(window, :) / trapz(binCenters, histograms(window, :));
    end
    
    disp("SSE umbrella sampling time below");
    toc;
    save(filenameWorkSpace);
    
    % Adjusting histogram weights using the TPS function and calculating SSE values
    [~, unbiased_histogram_SSE] = adjust_histogram_weightsTPS(normalized_histograms, binEdges, Overlaps);
    Ct_Umbrella_SSE = trapz(binCenters(BRange), unbiased_histogram_SSE(BRange));
    unbiased_free_energy_SSE = -kT * log(unbiased_histogram_SSE);
    unbiased_free_energy_SSE = unbiased_free_energy_SSE - min(unbiased_free_energy_SSE); % Shift to zero minimum free energy
    
    clear pathIn_window PsiIn_window;
    save(filenameWorkSpace);
    
    %% Unconstrained Sampling
    nPathsPerWindow = 3*nPathsPerWindow;
    filenameWorkSpace = sprintf('SymDeepworkspaceSSEUnconstrained_T=%.2f_Gamma=%.2f.mat', kT, Gamma);
    
    % Parameters
    tf = 40; % Total time duration
    trajectorySteps = round(tf / dt); % Number of SSE steps for new path
    TSpan = linspace(0, tf, trajectorySteps);
    nCopy = 16; % Number of parallel copies for sampling
    
    acceptRatio = zeros(1, nCopy);
    numGroups = 200; % Groups for statistical analysis
    totalPaths = zeros(1, nCopy);
    
    tic;
    h_AB_win = zeros(nCopy, numGroups);
    
    parfor copy = 1:nCopy
        count = 0;
        accepted = 0;
        path = zeros(1, trajectorySteps);
        Psi = zeros(bSize, trajectorySteps);
        stuck = 0;
    
        % Loop until path meets condition for regions A and B
        while path(1) > A(2) || ~any(path > B(1))
            eta = randn(1, trajectorySteps);
            PsiPert = DisplacementOpG(momentumPertUnc * ( 1i * randn) / sqrt(2), a, PsiGround);
            [path, Psi] = SSEDynamicsTPS(PsiPert, dt, tf, eta, HhatGamma, Lhat, Xhat, hbar, round(trajectorySteps / 2));
        end
    
        % Sampling paths and applying TPS shooting method
        while accepted < nPathsPerWindow + n_equilib
            PaccString = sprintf('Repeat');
            count = count + 1;
            stuck = stuck + 1;
            idx = ShootPointRND(path); % Choose random shooting point
            eta = randn(1, trajectorySteps);
            PsiPert = DisplacementOp(momentumPertUnc * ( 1i * randn) / sqrt(2), a) * Psi(:, idx);
            [newPath, newPsi] = SSEDynamicsTPS(PsiPert, dt, tf, eta, HhatGamma, Lhat, Xhat, hbar, idx);
            [choice, Ns] = PathCheckerUnconstrained(newPath, A, B);
%             plot(newPath) % uncomment to see path (also remove parfor)
            % Path acceptance and transformations
            if choice ~= 0
                [newPath, newPsi] = PathTransformerSSE(newPath, newPsi, choice);
                [revPath, ~] = PathTransformerSSE(path, Psi, choice);
                [~, barNs] = PathCheckerUnconstrained(revPath, A, B);
                Pac = PaccSSETransformed(newPsi, Psi, idx, ThermalState, HhatGamma, Lhat, hbar, dt, choice);
                Pac = (Ns / barNs) * Pac;
                PaccString = sprintf('Pacc_%s=%.2g', ChoiceString{choice}, Pac);

                % Accept path based on probability and reset "stuck" counter
                if rand < Pac
                    accepted = accepted + 1;
                    path = newPath;
                    Psi = newPsi;
                    stuck = 0;
                end
            end
    
            % Update histograms after equilibrium is reached
            if accepted > n_equilib && mod(count, tickrate) == 0
                totalPaths(copy) = totalPaths(copy) + 1;
                h_AB_win(copy, :) = h_AB_win(copy, :) + mean(reshape(indicatorB(path, B(2), B(1)), [], numGroups), 1);
            end
    
            % Display progress and current acceptance ratio
            fprintf('SSE copy %d, Progress: %5.2f, Stuck: %5.d, %s\n', copy , (accepted - n_equilib) / nPathsPerWindow, stuck, PaccString);
        end
    
        % Finalize path statistics for current copy
        h_AB_win(copy, :) = h_AB_win(copy, :) / totalPaths(copy);
        acceptRatio(copy) = accepted / totalPaths(copy);
        fprintf('SSE copy %d, finished with accept ratio: %.2f\n', copy, acceptRatio(copy));
    end
    
    fprintf('SSE unconstrained mean accept ratio: %.2g\n', mean(acceptRatio));
    
    disp("SSE unconstrained sampling time below");
    toc;
    save(filenameWorkSpace);
    
    %% Processing Results
    h_Star_AB_SSE = mean(h_AB_win, 1);
    dt_smooth = tf / numGroups;
    N_t_smooth = round(N_t_prime * (numGroups / trajectorySteps));
    h_Star_AB_SSE_Smooth = mean(reshape(h_Star_AB_SSE, [], numGroups), 1);
    
    % Calculating derivatives and correlation rate for SSE
    smoothSpan = linspace(0, tf, numGroups);
    dh_Star_AB_SSE = central_diff_8th_order(smoothSpan, h_Star_AB_SSE);
    Correl_SSE = h_Star_AB_SSE_Smooth / h_Star_AB_SSE_Smooth(N_t_smooth) * Ct_Umbrella_SSE;
    rate_SSE = dh_Star_AB_SSE / h_Star_AB_SSE_Smooth(N_t_smooth) * Ct_Umbrella_SSE;
    save(filenameWorkSpace);
end

%% Functions

% Displacement operator function for vector of displacements
function result = DisplacementOp(alphaVec, a)
    N = length(alphaVec);
    bSize = size(a, 1);
    result = zeros(bSize, bSize, N);
    for n = 1:N
        alpha = alphaVec(n);
        result(:, :, n) = expm(alpha * a' - conj(alpha) * a);
    end
end

% Displacement operator function with ground state
function result = DisplacementOpG(alphaVec, a, PsiGround)
    N = length(alphaVec);
    bSize = size(a, 1);
    result = zeros(bSize, N);
    for n = 1:N
        alpha = alphaVec(n);
        result(:, n) = expm(alpha * a' - conj(alpha) * a) * PsiGround;
    end
end
