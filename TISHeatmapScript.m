%% Data Entry & Parallel Setup
% license('test', 'Distrib_Computing_Toolbox')
% % Manually specify the number of cores requested
% numCores = 128; % Update this with the number of cores requested in your PBS script
% fprintf('Number of allocated cores: %d\n', numCores);
% 
% % Set up a parallel pool with the specified number of cores
% poolobj = gcp('nocreate');
% if ~isempty(poolobj)
%     delete(poolobj); % Delete existing pool if it exists
% end
% parpool('local', numCores); % Create a new parallel pool
% 
% % Display the number of workers in the parallel pool
% currentPool = gcp;
% fprintf('Number of workers in the current parallel pool: %d\n', currentPool.NumWorkers);

%% Hamiltonian Setup
% Define parameters for Hamiltonian and potential calculations
hbar = 1; % Reduced Planck's constant
Xview = 6;
h2 = 0.35; % Coefficient for the quadratic term in the potential
h4 = 0.01; % Coefficient for the quartic term in the potential
Xmin = -sqrt(h2 / (2 * h4)); % Minimum position in the potential
Pmin = 0;
Vact = -h4 * Xmin^4 + h2 * Xmin^2; % Potential energy at minimum

% Define well boundaries and time parameters
A = [-6; -4];
B = flip(abs(A)); % Right well boundaries
tf = 200; % Final time for simulations

%% Transition Rate Setup
% Define parameters for path sampling and averaging
nPathsPerInterface = 2400; % Number of paths per interface
n_equilib = 200; % Number of equilibration paths
nAvg = 10; % Number of averages per temperature/gamma pair

% Temperature and dissipation ranges for the simulations
TempRange = Vact .* (0.1:0.04:0.5);
nSweep = length(TempRange); % Number of temperature/gamma values to sweep
GammaRange = linspace(0.1, 0.5, nSweep); % Range of Gamma values
totalIterations = nSweep^2 * nAvg; % Total number of iterations to run
tic % Start timing the execution

% Initialize rates
RateLangevinLinear=nan(1,nSweep^2*nAvg);
RateSSELinear=nan(1,nSweep^2*nAvg);

% % Load existing output if continuing from a previous run
% filename = sprintf('Robson_SingleTIS_job.o9187478');
% [RateLangevinLinear, RateSSELinear] = ReadOutputIncomplete(filename, TempRange, GammaRange, nAvg);

%% Parallelized Transition Rate Calculation
% Main parallel loop to calculate transition rates
parfor idx = 1:totalIterations
    % Determine the averaging index and temperature/gamma indices
    avgIdx = mod(idx-1, nAvg) + 1; % Averaging iteration
    trueIdx = ceil(idx / nAvg); % Overall index for temperature/gamma pairs
    m = mod(trueIdx-1, nSweep) + 1; % Temperature index
    k = ceil(trueIdx / nSweep); % Gamma index
    
    % Fetch current transition rates if they exist
    RateSSE = RateSSELinear(idx);
    RateLangevin = RateLangevinLinear(idx);
    
    % Call the function to compute rates if missing, or use existing values
    [RateSSELinear(idx), RateLangevinLinear(idx)] = computeRates( ...
        idx, nAvg, nSweep, GammaRange, TempRange, h4, h2, ...
        nPathsPerInterface, n_equilib, tf, A, RateLangevin, RateSSE);
end

% Save the workspace after the sweep is completed
save('workspaceTransitionRateTISHeatmapFinal.mat')
disp('Real Sweep time below')
toc % Stop timing the execution

%% Organize Data for Averaging
% Reshape data for averaging across different trials
RateSSECombined = reshape(RateSSELinear, [nAvg, totalIterations/nAvg]);
RateLangevinCombined = reshape(RateLangevinLinear, [nAvg, totalIterations/nAvg]);

% Compute the average rates across trials
averageRateSSECombined = reshape(mean(RateSSECombined, 1), [nSweep, nSweep]);
averageRateLangevinCombined = reshape(mean(RateLangevinCombined, 1), [nSweep, nSweep]);
save('workspaceTransitionRateTISHeatmapFinal.mat') % Save averaged data

% Close the parallel pool
delete(gcp);

%% Helper Function for Rate Calculation
function [RateSSEout, RateLangevinout] = computeRates(idx, nAvg, nSweep, GammaRange, TempRange, h4, h2, nPathsPerInterface, n_equilib, tf, A, RateLangevin, RateSSE)
    % Computes SSE and Langevin transition rates for given parameters
    % Inputs:
    % - idx: Current index in the loop
    % - nAvg, nSweep: Averaging count and number of temperature/gamma values
    % - GammaRange, TempRange: Ranges for dissipation and temperature
    % - h4, h2: Hamiltonian parameters
    % - nPathsPerInterface, n_equilib: Path sampling parameters
    % - tf: Final time for simulations
    % - A: Left well boundary positions
    % - RateLangevin, RateSSE: Initial transition rates (if any)
    % Outputs:
    % - RateSSEout, RateLangevinout: Computed or reused transition rates

    % Calculate indices for temperature and dissipation constant
    trueIdx = ceil(idx / nAvg); % Overall index
    m = mod(trueIdx-1, nSweep) + 1; % Temperature index
    k = ceil(trueIdx / nSweep); % Gamma index
    
    % Set dissipation and temperature for the current run
    Gamma = GammaRange(k);
    kT = TempRange(m);
    
    % Compute Langevin rate if not previously calculated
    if isnan(RateLangevin)
        RateLangevinout = LangevinTransitionRateTIS(kT, Gamma, h4, h2, nPathsPerInterface, n_equilib, tf, A);
    else
        RateLangevinout = RateLangevin;
    end
    fprintf('T_rate_Langevin: %.4g,   kT: %.4g,   Gamma: %.4g\n', RateLangevinout, kT, Gamma);
    
    % Compute SSE rate if not previously calculated
    if isnan(RateSSE)
        RateSSEout = SSETransitionRateTIS(kT, Gamma, h4, h2, nPathsPerInterface, n_equilib, tf, A);
    else
        RateSSEout = RateSSE;
    end
    fprintf('T_rate_SSE: %.4g,   kT: %.4g,   Gamma: %.4g\n', RateSSEout, kT, Gamma);
end
