function T_rate=SSETransitionRateTIS(kT,Gamma,h4,h2,nPathsPerInterface,n_equilib,tf,A)
% Calculates the transition rate using Stochastic Schrödinger Equation (SSE) dynamics and Transition Interface Sampling (TIS).
% Inputs:
% kT: Temperature.
% Gamma: Dissipation constant.
% h4, h2: Coefficients of the quartic and quadratic terms in the potential energy.
% nPathsPerInterface: Number of unique paths per interface.
% n_equilib: Number of equilibration paths.
% tf: Final time for path evolution.
% A: Interface positions for the left well.

Xmin=-sqrt(h2/2/h4); % Minimum position in the potential well.
B=flip(abs(A)); % Flips and takes absolute value of A for right well interface positions.
hbar=1; % Planck's constant over 2*pi, set to 1 for simplicity.
bSize=30; % Basis size for the quantum harmonic oscillator.
a=zeros(bSize); % Annihilation operator in matrix form.
for n=1:bSize-1
    a(n,n+1)=sqrt(n); % Constructing the annihilation operator.
end

Xhat=sqrt(hbar/2)*(a'+a); % Position operator.
Phat=1i*sqrt(hbar/2)*(a'-a); % Momentum operator.
Lhat=sqrt(4*Gamma*kT/hbar)*Xhat+1i*sqrt(Gamma*hbar/4/kT)*Phat; % Coupling to thermal bath.

HhatGamma=0.5*Phat^2+h4*Xhat^4-h2*Xhat^2+0.5*Gamma*(Xhat*Phat+Phat*Xhat); % Hamiltonian including the system-bath interaction.
ThermalState=LindbladNull(HhatGamma,Lhat,hbar); % Calculates the thermal steady state.
PsiGround=zeros(bSize,1); % Ground state initialization.
PsiGround(1)=1; % Ground state is the first basis vector.
dt=0.001;
tickrate=3; % Sampling frequency for paths.

%% T loop minus - Sampling in the left well
% This section initializes path variables and performs SSE dynamics to sample paths starting from the left well, targeting transitions to the right.

clear path Psi % Clears previous path and wavefunction information.
LeftWell= DisplacementOpG(Xmin/sqrt(2),a,PsiGround); % Initial state preparation in the left well.

% Loop to generate and accept paths based on the transition probability.
accepted = 0; % Counter for accepted paths.
check=0; % Flag to check if a valid path has been generated.
path=nan; % Initializes path variable.
Psi=nan; % Initializes wavefunction variable.
Pac=nan; % Initializes acceptance probability variable.
momentumPert=0.3*rand; % Random momentum perturbation for path variation.
while check==0 % Loop to ensure at least one valid path is generated.
    [path, Psi, ~,check,~] =SSEDynamicsTISLeft(LeftWell,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2)); % Generates a path starting from the left well.
end
% Initializes variables for tracking total time and paths.
TotalTime=0;
TotalPaths=0;
ticker=0; % Counter for thinning.
stuck=0; % Counter for stuck paths.

while accepted -n_equilib< nPathsPerInterface % Loop to collect the required number of paths.
    ticker=ticker+1;
    stuck=stuck+1;
    idx= ShootPointRND(path); % Chooses a random shooting point along the path.
    PsiPert=DisplacementOp(1i*momentumPert*randn/sqrt(2),a)*Psi(:,idx); % Applies a momentum perturbation.
    [newPath,newPsi,idx_new,check,~] =SSEDynamicsTISLeft(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2)); % Generates a new path from the perturbed state.
    if check==1 % Checks if the new path is valid.
        Pac=PaccSSETIS(newPsi,Psi,idx,idx_new,ThermalState,HhatGamma,Lhat,hbar,dt); % Calculates the acceptance probability.
        if rand <Pac % Accepts or rejects the new path based on Pac.
            stuck=0;
            accepted = accepted + 1;
            Psi=newPsi;
            path=newPath;
        end
    end
    if accepted>n_equilib&&mod(ticker,tickrate)==0 % Thinning to reduce autocorrelation.
        TotalTime=TotalTime+length(path)*dt;
        TotalPaths=TotalPaths+1;
    end
    % Display progress information.
    % fprintf('SSE Left: Progress: %5.2f,  Steps: %15.3g,   Pac: %15.2g, TAvg:%15.2g, Stuck: %5.d\n',(accepted-n_equilib)/nPathsPerInterface, length(path), Pac,TotalTime/TotalPaths,stuck)
end
AvgTMinus=TotalTime/TotalPaths; % Average time for paths in the left well.

%% Sampling SSE Right TIS Variable Windows
% This section samples paths in the to the right of the first interface using variable window sizes to ensure efficient sampling across the transition space.

Prob=0.3+0.3*rand; % Probability for determining next window size.
clear XVals path Psi % Clears variables for right well sampling.
XVals=A(2); % Initial interface position for the right well.
binPositions = A(2):0.025:B(1); % Defines bin positions for histogram.
nBins = length(binPositions); % Number of bins.
CrossingHistograms=zeros(1,nBins); % Initializes crossing histograms.
successfrac=0; % Success fraction initialization.
TimeTracker=0; % Time tracker initialization.
converged=false; % Flag for convergence.
window=1; % Current window index.

while converged==false % Loop until convergence criteria are met.
    accepted = 0; % Counter for accepted paths.
    TotalTime=0; % Total time tracker.
    check=0; % Flag for valid path generation.
    path=nan; % Initializes path variable.
    Psi=nan; % Initializes wavefunction variable.
    Pac=nan; % Initializes acceptance probability variable.
    stuck=0; % Counter for stuck paths.
    TempHist=zeros(1,nBins); % Temporary histogram for current window.
    TempMP=[]; % Temporary max positions tracker.
    ticker=0; % Counter for thinning.

    while check==0 % Loop to generate at least one valid path for the current window.
        PsiPert= DisplacementOpG((XVals(window)+0.1+1i*momentumPert*randn)/sqrt(2),a,PsiGround); % Perturbs initial state.
        [path, Psi, ~,check,~] =SSEDynamicsTISRight(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2),B(1)); % Generates a path.
    end

    while (accepted-n_equilib)< nPathsPerInterface % Loop to collect the required number of paths for current window.
        ticker=ticker+1;
        stuck=stuck+1;
        idx= ShootPointRND(path); % Chooses a shooting point based on interface position.
        PsiPert=DisplacementOp(1i*momentumPert*randn/sqrt(2),a)*Psi(:,idx); % Applies momentum perturbation.
        [newPath, newPsi, idx_new,check,~] =SSEDynamicsTISRight(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2),B(1)); % Generates a new path.
        MP=max(newPath); % Finds the maximum position in the new path.
        if check==1&&MP>=XVals(window) % Checks if the new path is valid and crosses the current window.
            Pac=PaccSSETIS(newPsi,Psi,idx,idx_new,ThermalState,HhatGamma,Lhat,hbar,dt); % Calculates the acceptance probability.
            if rand <Pac % Accepts or rejects the new path based on Pac.
                path=newPath;
                Psi=newPsi;
                accepted=accepted+1;
                stuck=0;
            end
        end
        MP=max(path); % Updates the maximum position for progress tracking.
        % Display progress information for the current window.
        % fprintf('Right Window: %5.1g, Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,   Pac: %15.2g, Stuck: %5.d\n', window,(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
        if accepted>n_equilib&&mod(ticker,tickrate)==0 % Thinning to reduce autocorrelation.
            TempMP=[TempMP, MP]; % Updates temporary max positions tracker.
            for n=1:length(MP) % Updates the histogram based on crossing events.
                bump = customHeaviside(binPositions - XVals(window)) - customHeaviside(binPositions - MP(n));
                TempHist=TempHist+bump;
            end
            TotalTime=TotalTime+dt*length(path); % Updates total time for the current window.
        end
    end
    CrossingHistograms(window,:)=TempHist/length(TempMP); % Updates crossing histograms.
    TimeTracker(window)=TotalTime/length(TempMP); % Updates time tracker.
    TempMP=sort(TempMP); % Sorts temporary max positions.
    XVals(window+1)=max(roundToNearest(TempMP(round((1-Prob)*length(TempMP))),0.025),XVals(window)+0.025); % Determines the next window position.
    if XVals(window+1)>=0.5 % Checks for convergence.
        converged=true;
        XVals(window+1)=B(1); % Sets final window position to the right well boundary.
        successfrac(window)=length(find(TempMP>=XVals(window+1)))/length(TempMP); % Calculates success fraction for the current window.
    else % If not converged, prepares for the next window.
        CrossingHistograms(window+1,:)=zeros(nBins,1); % Initializes histogram for the next window.
        TimeTracker=[TimeTracker, 0]; % Extends time tracker for the next window.
        successfrac(window)=length(find(TempMP>=XVals(window+1)))/length(TempMP); % Updates success fraction for the current window.
        % Display information for the completed window and setup for the next.
        window=window+1; % Advances to the next window.
    end
end
clear TempMP TempHist % Clears temporary variables.
[adjusted_weights,~] =adjust_histogram_weightsTIS(CrossingHistograms,XVals,binPositions,successfrac); % Adjusts histogram weights.

% Calculates average time for paths crossing from left to right.
AvgTPlus=sum(adjusted_weights.'.*TimeTracker);
AvgTLoop=AvgTPlus+AvgTMinus; % Total average time for complete transition.
PTran=prod(successfrac); % Product of success fractions.
NTries=1/PTran; % Number of tries to achieve a successful transition.
T_rate=1/(NTries*AvgTLoop); % Transition rate calculation.

end

%% Supporting functions
% These functions provide utility operations such as state displacement, custom Heaviside function, and rounding to nearest value.

function result = DisplacementOp(alphaVec,a)
% Generates displacement operators for a given set of complex amplitudes.
% Inputs:
% alphaVec: Vector of complex amplitudes for the displacement.
% a: Annihilation operator matrix.
% Output:
% result: Tensor of displacement operators for each amplitude in alphaVec.

N = length(alphaVec);
bSize=size(a,1);
result = zeros(bSize, bSize, N);
for n = 1:N
    alpha = alphaVec(n);
    result(:,:,n) = expm(alpha*a' - conj(alpha)*a); % Exponential of the displacement operator.
end
end

function result = DisplacementOpG(alphaVec,a,PsiGround)
% Applies displacement operators to the ground state.
% Inputs:
% alphaVec: Vector of complex amplitudes for the displacement.
% a: Annihilation operator matrix.
% PsiGround: Ground state vector.
% Output:
% result: Matrix of displaced ground states for each amplitude in alphaVec.

N = length(alphaVec);
bSize=size(a,1);
result = zeros(bSize, N);
for n = 1:N
    alpha = alphaVec(n);
    result(:,n) = expm(alpha*a' - conj(alpha)*a)*PsiGround; % Displaced ground state.
end
end

function y = customHeaviside(x)
% Custom Heaviside step function.
% Inputs:
% x: Input value or array.
% Output:
% y: Heaviside step function applied to x (1 for x >= 0, 0 otherwise).

y = double(x >= 0);  % Returns 1 for x >= 0 and 0 otherwise.
end

function roundedValue = roundToNearest(value,num)
% Rounds a value to the nearest multiple of a given number.
% Inputs:
% value: The value to be rounded.
% num: The number to which value will be rounded to the nearest multiple.
% Output:
% roundedValue: The rounded value.

roundedValue = round(value / num) * num; % Rounded to the nearest multiple of num.
end
