function T_rate=LangevinTransitionRateTIS(kT,Gamma,h4,h2,nPathsPerInterface,n_equilib,tf,A)
% Calculates the transition rate for a system under Langevin dynamics using the Transition Interface Sampling (TIS) approach.
% Inputs:
% kT: Temperature.
% Gamma: Dissipation constant.
% h4, h2: Coefficients of the quartic and quadratic terms in the potential energy.
% nPathsPerInterface: Number of unique paths to sample per interface.
% n_equilib: Number of equilibration paths.
% tf: Final simulation time.
% A: Array of interface positions for path sampling in the left well.

Xmin=-sqrt(h2/2/h4); % Minimum of the double-well potential.
B=flip(abs(A)); % Interface positions for path sampling in the right well.
dt=0.001; % Time step for the Langevin dynamics simulation.
tickrate=3; % Sampling frequency for paths.

%% T loop minus - Sampling paths starting in the left well
tic % Starts timer for performance measurement.
accepted = 0; % Counter for accepted paths.
check=0; % Flag to check if a valid path has been generated.
TotalPaths=0; % Counter for total paths sampled.
ticker=0; % Counter for sampling rate control.
path=nan; % Initializes the path variable.
Psi=nan; % Initializes the state vector variable.
Pac=nan; % Initializes the acceptance probability variable.
stuck=0; % Counter for the number of iterations without accepting a new path.
momentumPert=0.3*rand; % Initial random momentum perturbation.
while check==0 % Initial loop to generate at least one valid path.
    [path, Psi, ~,check,~] = LangevinDynamicsTISLeft([Xmin; sqrt(kT)*randn],dt,tf,Gamma,kT,h4,h2,A(2));
end
TotalTime=0; % Initializes total simulation time for accepted paths.

while accepted -n_equilib< nPathsPerInterface % Sampling loop for the left well.
    ticker=ticker+1;
    stuck=stuck+1;
    idx= ShootPointRND(path); % Randomly chooses a shooting point along the path.
    PsiPert=Psi(:,idx); % Extracts the state at the chosen shooting point.
    PsiPert(2)=PsiPert(2)+momentumPert*randn; % Applies a momentum perturbation.
    [newPath,newPsi,idx_new,check,~] =LangevinDynamicsTISLeft(PsiPert,dt,tf,Gamma,kT,h4,h2,A(2));
    if check==1 % Checks if the new path is valid.
        Pac=PaccLangevinTIS(newPsi,Psi,idx,idx_new,kT,Gamma,dt,h4,h2); % Calculates the acceptance probability.
        if rand <Pac % Accepts or rejects the new path based on Pac.
            accepted = accepted + 1;
            Psi=newPsi;
            path=newPath;
            stuck=0;
        end
    end
    if accepted>n_equilib&&mod(ticker,tickrate)==0 % Calculates total time and paths for thinning.
        TotalTime=TotalTime+length(path)*dt;
        TotalPaths=TotalPaths+1;
    end
    % Displays progress information for sampling in the left well.
    % fprintf('Langevin Left: Progress: %5.2f,  Steps: %15.3g,   Pac: %15.2g, TAvg:%15.2g, Stuck: %5.d\n',(accepted-n_equilib)/nPathsPerInterface, length(path), Pac,TotalTime/TotalPaths,stuck)
end
AvgTMinus=TotalTime/TotalPaths; % Average time for paths starting in the left well.

%% Initial Right TIS Variable Windows
% This section samples paths starting in the right well using variable windows to ensure efficient sampling.

Prob=0.3+0.3*rand; % Random probability for determining next window size.
clear XVals path Psi % Clears variables for right well sampling.
XVals=A(2); % Initial interface position for the right well sampling.
binPositions = A(2):0.025:B(1); % Defines bin positions for histogram calculation.
nBins = length(binPositions); % Number of bins for histogram.
CrossingHistograms=zeros(1,nBins); % Initializes crossing histograms.
successfrac=0; % Initializes success fraction variable.
TimeTracker=0; % Initializes time tracker variable.
converged=false; % Flag for convergence in the variable window sampling.
window=1; % Current window index for right well sampling.

while converged==false % Loop until convergence criteria are met.
    accepted = 0; % Counter for accepted paths.
    TotalTime=0; % Total time tracker for current window.
    check=0; % Flag for valid path generation.
    ticker=0; % Counter for thinning control.
    path=nan; % Initializes path variable for right well sampling.
    Psi=nan; % Initializes state vector variable for right well sampling.
    stuck=0; % Counter for the number of iterations without accepting a new path in the right well.
    Pac=nan; % Initializes the acceptance probability variable for right well sampling.
    while check==0 % Initial loop to generate at least one valid path for the current window in the right well.
        [path, Psi, ~,check,~] = LangevinDynamicsTISRight([XVals(window)+0.01; sqrt(kT)*randn],dt,tf,Gamma,kT,h4,h2,A(2),B(1));
    end
    TempHist=zeros(1,nBins); % Temporary histogram for current window.
    TempMP=[]; % Temporary max positions tracker for current window.
    while (accepted-n_equilib)< nPathsPerInterface % Sampling loop for the current window in the right well.
        ticker=ticker+1;
        stuck=stuck+1;
        idx= ShootPointRND(path); % Chooses a shooting point based on interface position.
        PsiPert=Psi(:,idx); % Extracts the state at the chosen shooting point.
        PsiPert(2)=PsiPert(2)+momentumPert*randn; % Applies a momentum perturbation.
        [newPath, newPsi, idx_new,check,~] = LangevinDynamicsTISRight(PsiPert,dt,tf,Gamma,kT,h4,h2,A(2),B(1));
        MP=max(newPath); % Finds the maximum position in the new path.
        if check==1&&MP>=XVals(window) % Checks if the new path is valid and crosses the current window.
            Pac=PaccLangevinTIS(newPsi,Psi,idx,idx_new,kT,Gamma,dt,h4,h2); % Calculates the acceptance probability.
            if rand <Pac % Accepts or rejects the new path based on Pac.
                path=newPath;
                Psi=newPsi;
                accepted=accepted+1;
                stuck=0;
            end
        end
        MP=max(path); % Updates the maximum position for progress tracking.
        % Displays progress information for the current window in the right well.
        % fprintf('Right Window: %5.1g, Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,   Pac: %15.2g, Stuck: %5.d\n', window,(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
        if accepted>n_equilib&&mod(ticker,tickrate)==0 % Thinning to reduce autocorrelation in the right well.
            TempMP=[TempMP, MP]; % Updates temporary max positions tracker for the current window.
            for n=1:length(MP) % Updates the histogram based on crossing events for the current window.
                bump = customHeaviside(binPositions - XVals(window)) - customHeaviside(binPositions - MP(n));
                TempHist=TempHist+bump;
            end
            TotalTime=TotalTime+dt*length(path); % Updates total time for the current window.
        end
    end
    CrossingHistograms(window,:)=TempHist/length(TempMP); % Updates crossing histograms for the current window.
    TimeTracker(window)=TotalTime/length(TempMP); % Updates time tracker for the current window.
    TempMP=sort(TempMP(:)); % Sorts temporary max positions for the current window.
    XVals(window+1)=roundToNearest(TempMP(round((1-Prob)*length(TempMP))),0.025); % Determines the next window position for the right well.
    if XVals(window+1)>=0.5 % Checks for convergence in the variable window sampling.
        converged=true;
        XVals(window+1)=B(1); % Sets final window position to the right well boundary.
        successfrac(window)=length(find(TempMP>=XVals(window+1)))/length(TempMP); % Calculates success fraction for the current window.
    else % If not converged, prepares for the next window in the right well.
        CrossingHistograms(window+1,:)=zeros(nBins,1); % Initializes histogram for the next window.
        TimeTracker=[TimeTracker, 0]; % Extends time tracker for the next window.
        successfrac(window)=length(find(TempMP>=XVals(window+1)))/length(TempMP); % Updates success fraction for the current window.
        % Display information for the completed window and setup for the next in the right well.
        window=window+1; % Advances to the next window in the right well.
    end
end
clear TempMP TempMP TempHist % Clears temporary variables for right well sampling.
[adjusted_weights,~] =adjust_histogram_weightsTIS(CrossingHistograms,XVals,binPositions,successfrac); % Adjusts histogram weights for the right well.
AvgTPlus=sum(adjusted_weights.'.*TimeTracker); % Calculates average time for paths starting in the right well.
AvgTLoop=AvgTPlus+AvgTMinus; % Total average time for complete transition.
PTran=prod(successfrac); % Product of success fractions for the transition.
NTries=1/PTran; % Number of tries to achieve a successful transition.
T_rate=1/(NTries*AvgTLoop); % Transition rate calculation for the system under Langevin dynamics with TIS.
end

% Additional functions for utility operations such as custom Heaviside function and rounding to nearest value

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

