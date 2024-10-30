%% Data entry
clearvars
nCopy = 16; % Update this with the number of cores requested in your PBS script

nPathsPerInterface=1800;
n_equilib=100;
tickrate=10;

% nPathsPerInterface=12;
% n_equilib=2;
%% Hamiltonian

Xview=8;
h2=0.35; %Quartic
h4=0.01;
Xmin=-sqrt(h2/2/h4);
Vact=-h4*Xmin^4+h2*Xmin^2;
Gamma=0.25;
kT=0.25*Vact;
dt=0.001;
tf=80;
VarMomentum=0.2;
VarTrand=0.05;
Prob=0.2;


%Quartic
hbar=1;

V = @(x) (h4*x.^4-h2*x.^2);
Vscaled = @(x) (h4*x.^4-h2*x.^2-V(Xmin)-kT);

H = @(Z) (0.5*Z(2)^2+h4*Z(1)^4-h2*Z(1)^2);
A=[-6;-4];

B=flip(abs(A));
% B=[3; 5];
% XVals=[-3 -1 3];


bSize=30; %basis size
a=zeros(bSize);
for n=1:bSize-1
    a(n,n+1)=sqrt(n);
end

Xhat=sqrt(hbar/2)*(a'+a);
Phat=1i*sqrt(hbar/2)*(a'-a);
Lhat=sqrt(4*Gamma*kT/hbar)*Xhat+1i*sqrt(Gamma*hbar/4/kT)*Phat; %caldeira


HhatOrig=0.5*Phat^2+h4*Xhat^4-h2*Xhat^2;
HhatGamma=0.5*Phat^2+h4*Xhat^4-h2*Xhat^2+0.5*Gamma*(Xhat*Phat+Phat*Xhat);
ThermalState=LindbladNull(HhatGamma,Lhat,hbar);
% ThermalState=expm((-1/kT)*HhatOrig);
% ThermalState=ThermalState/trace(ThermalState);
PsiGround=zeros(bSize,1); %Initial state
PsiGround(1)=1;
filenameWorkSpace = sprintf('Cl2usterWorkspaceTestFluxAndProbabilities_Two_T=%.2f_Gamma=%.2f.mat', kT,Gamma);


%% Langevin Left
tic
% Prob=ProbIdeal;
clear XValsLangevinLeft path Psi
XValsLangevinLeft=A(2);
binPositionsLeft = A(2):-0.025:-8;
nBinsLeft = length(binPositionsLeft);
CrossingHistogramsLangevinLeft=zeros(nBinsLeft,nCopy,1);
successfracLangevinLeft=0;
TimeTrackerLangevinLeft=zeros(nCopy,1);
TransitionTimeLangevinLeft=zeros(nCopy,1);
converged=false;
window=1;
% Loop LangevinLeft
while converged==false
    MPArrayTemp=cell(1,nCopy);
    parfor copy = 1:nCopy
        rng(100*copy + sum(clock)*1000);  % Set unique seed based on clock and loop index
        accepted = 0;
        success=0;
        TotalTime=0;
        check=0;
        ticker=0;
        bump=zeros(1,nBinsLeft);
        path=nan;
        Psi=nan;
        MP=nan;
        stuck=0;
        while check==0
            [path, Psi, ~,check,~] = LangevinDynamicsTISLeft([XValsLangevinLeft(window)-0.01; sqrt(kT)*randn],dt,tf,Gamma,kT,h4,h2,A(2));
        end
        TempHist=zeros(1,nBinsLeft);
        TempMP=[];
        while (accepted-n_equilib)< nPathsPerInterface
            Pac=-1;
            ticker=ticker+1;
            stuck=stuck+1;
            idx= ShootPointRND(path); % Choose shooting point
            PsiPert=Psi(:,idx);
            PsiPert(2)=PsiPert(2)+VarMomentum*randn;
            [newPath, newPsi, idx_new,check,~] = LangevinDynamicsTISLeft(PsiPert,dt,tf,Gamma,kT,h4,h2,A(2));
            MP=min(newPath);
            if check==1&&MP<=XValsLangevinLeft(window)
                Pac=PaccLangevinTIS(newPsi,Psi,idx,idx_new,kT,Gamma,dt,h4,h2);
                if rand <Pac
                    path=newPath;
                    Psi=newPsi;
                    accepted=accepted+1;
                    stuck=0;
                end
            end
            MP=min(path);
            fprintf('Langevin Left Window: %5.d, Position %5.4g,  Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,   Pac: %15.2g, Stuck: %5.d\n', window,XValsLangevinLeft(window),(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
            if accepted>n_equilib&&mod(ticker,tickrate)==0
                TempMP=[TempMP, MP];
                for n=1:length(MP)
                    bump = abs(customHeavisideLeft(binPositionsLeft - XValsLangevinLeft(window)) - customHeavisideLeft(binPositionsLeft - MP(n)));
                    TempHist=TempHist+bump;
                end
                TotalTime=TotalTime+dt*length(path);
            end
        end
        MPArrayTemp{copy}=TempMP;
        CrossingHistogramsLangevinLeft(:,copy,window)=TempHist/length(TempMP);
        TimeTrackerLangevinLeft(copy,window)=TotalTime/length(TempMP);
    end
    MPArrayTemp=horzcat(MPArrayTemp{:});
    MPArrayTemp=sort(MPArrayTemp(:),"descend");
    MPArrayLangevinLeftWindow{window}=MPArrayTemp;
    XValsLangevinLeft(window+1)=min(roundToNearest(MPArrayTemp(round((1-Prob)*length(MPArrayTemp))),0.025),XValsLangevinLeft(window)-0.025); % Determines the next window position.
    successfracLangevinLeft(window)=length(find(MPArrayTemp<=XValsLangevinLeft(window+1)))/length(MPArrayTemp);
    if window==2
        converged=true;
    else
        CrossingHistogramsLangevinLeft(:,:,window+1)=zeros(nBinsLeft,nCopy);
        TimeTrackerLangevinLeft=[TimeTrackerLangevinLeft, zeros(nCopy,1)];
        fprintf('Window: %5.1g finished,     Next window: %5.2f,    successfrac: %15.3g\n', window,XValsLangevinLeft(window+1), successfracLangevinLeft(window))
        window=window+1;
    end
end
fprintf('Window: %5.1g finished, successfrac: %15.3g\n', window, successfracLangevinLeft(window))
clear MPArrayTemp TempMP TempHist
CrossingHistogramsLangevinLeft=squeeze(mean(CrossingHistogramsLangevinLeft,2)).';
[adjusted_weightsLangevinLeft,unbiasedhistogramLangevinLeft] =adjust_histogram_weightsTIS(CrossingHistogramsLangevinLeft,abs(XValsLangevinLeft),abs(binPositionsLeft),successfracLangevinLeft);
AvgTMinusLangevin=sum(repmat(adjusted_weightsLangevinLeft.',nCopy,1).*TimeTrackerLangevinLeft,2);
% AvgTMinusLangevin=TimeTrackerLangevinLeft(:,1);


toc

%% Langevin Right

tic
% Prob=ProbIdeal;
clear XValsLangevin path Psi
XValsLangevin=A(2);
binPositions = A(2):0.025:B(1);
nBins = length(binPositions);
CrossingHistogramsLangevin=zeros(nBins,nCopy,1);
successfracLangevin=0;
TimeTrackerLangevin=zeros(nCopy,1);
TransitionTimeLangevin=zeros(nCopy,1);
converged=false;
window=1;
% Loop Langevin
while converged==false
    MPArrayTemp=cell(1,nCopy);
    parfor copy = 1:nCopy
        rng(100*copy + sum(clock)*1000);  % Set unique seed based on clock and loop index
        accepted = 0;
        success=0;
        TotalTime=0;
        check=0;
        ticker=0;
        bump=zeros(1,nBins);
        path=nan;
        Psi=nan;
        MP=nan;
        stuck=0;
        while check==0
            [path, Psi, ~,check,~,~] = LangevinDynamicsTISRight([XValsLangevin(window)+0.01; sqrt(kT)*randn],dt,tf,Gamma,kT,h4,h2,A(2),B(1));
        end
        TempHist=zeros(1,nBins);
        TempMP=[];
        while (accepted-n_equilib)< nPathsPerInterface
            Pac=-1;
            ticker=ticker+1
            stuck=stuck+1;
            idx= ShootPointRND(path); % Choose shooting point
            PsiPert=Psi(:,idx);
            PsiPert(2)=PsiPert(2)+VarMomentum*randn;
            [newPath, newPsi, idx_new,check,~,~] = LangevinDynamicsTISRight(PsiPert,dt,tf,Gamma,kT,h4,h2,A(2),B(1));
            MP=max(newPath);
            if check==1&&MP>=XValsLangevin(window)
                Pac=PaccLangevinTIS(newPsi,Psi,idx,idx_new,kT,Gamma,dt,h4,h2);
                if rand <Pac
                    path=newPath;
                    Psi=newPsi;
                    accepted=accepted+1;
                    stuck=0;
                end
            end
            MP=max(path);
            fprintf('Langevin Right Window: %5.d, Position %5.4g,  Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,   Pac: %15.2g, Stuck: %5.d\n', window,XValsLangevin(window),(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
            if accepted>n_equilib&&mod(ticker,tickrate)==0
                TempMP=[TempMP, MP];
                for n=1:length(MP)
                    bump = customHeaviside(binPositions - XValsLangevin(window)) - customHeaviside(binPositions - MP(n));
                    TempHist=TempHist+bump;
                end
                TotalTime=TotalTime+dt*length(path);
            end
        end
        MPArrayTemp{copy}=TempMP;
        CrossingHistogramsLangevin(:,copy,window)=TempHist/length(TempMP);
        TimeTrackerLangevin(copy,window)=TotalTime/length(TempMP);
    end
    MPArrayTemp=horzcat(MPArrayTemp{:});
    MPArrayTemp=sort(MPArrayTemp(:));
    MPArrayLangevinWindow{window}=MPArrayTemp;
    XValsLangevin(window+1)=max(roundToNearest(MPArrayTemp(round((1-Prob)*length(MPArrayTemp))),0.025),XValsLangevin(window)+0.025); % Determines the next window position.
    if XValsLangevin(window+1)>=0.5
        converged=true;
        XValsLangevin(window+1)=B(1);
        successfracLangevin(window)=length(find(MPArrayTemp>=XValsLangevin(window+1)))/length(MPArrayTemp);
    else
        CrossingHistogramsLangevin(:,:,window+1)=zeros(nBins,nCopy);
        TimeTrackerLangevin=[TimeTrackerLangevin, zeros(nCopy,1)];
        successfracLangevin(window)=length(find(MPArrayTemp>=XValsLangevin(window+1)))/length(MPArrayTemp);
        fprintf('Window: %5.1g finished,     Next window: %5.2f,    successfrac: %15.3g\n', window,XValsLangevin(window+1), successfracLangevin(window))
        window=window+1;
    end
end
fprintf('Window: %5.1g finished, successfrac: %15.3g\n', window, successfracLangevin(window))
clear MPArrayTemp TempMP TempHist
CrossingHistogramsLangevin=squeeze(mean(CrossingHistogramsLangevin,2)).';
[adjusted_weightsLangevin,unbiasedhistogramLangevin] =adjust_histogram_weightsTIS(CrossingHistogramsLangevin,XValsLangevin,binPositions,successfracLangevin);
AvgTPlusLangevin=sum(repmat(adjusted_weightsLangevin.',nCopy,1).*TimeTrackerLangevin,2);
% AvgTPlusLangevin=TimeTrackerLangevin(:,1);
AvgTLoopHistLangevin=AvgTPlusLangevin+AvgTMinusLangevin;

PTranLangevin=prod(successfracLangevin);
fprintf('PTranLangevin %15.3g\n', PTranLangevin)

toc

%% Langevin rates
clear  pathInLangevin PsiInLangevin  pathInLangevinLeft PsiInLangevinLeft Psi path

N_TriesLangevin=1./prod(successfracLangevin,2);
T_rateLangevinHistLoop=1/mean((N_TriesLangevin).*AvgTLoopHistLangevin)

%% SSE Left
tic
clear XValsSSELeft path Psi
XValsSSELeft=A(2);
window=1;
binPositionsLeft = A(2):-0.0125:-7;
nBinsLeft = length(binPositionsLeft);
CrossingHistogramsSSELeft=zeros(nBinsLeft,nCopy,1);
successfracSSELeft=0;
TimeTrackerSSELeft=zeros(nCopy,1);
converged=false;
%% Loop
while converged==false
    MPArrayTemp=cell(1,nCopy);
    parfor copy = 1:nCopy
        rng(100*copy + sum(clock)*1000);  % Set unique seed based on clock and loop index
        accepted = 0;
        success=0;
        TotalTime=0;
        check=0;
        path=nan;
        Psi=nan;
        stuck=0;
        TempHist=zeros(1,nBinsLeft);
        TempMP=[];
        ticker=0;
        while check==0
            PsiPert= DisplacementOpG((XValsSSELeft(window)-0.2+1i*VarMomentum*randn)/sqrt(2),a,PsiGround);
            [path, Psi, ~,check,~] =SSEDynamicsTISLeft(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2));
        end
        while (accepted-n_equilib)< nPathsPerInterface
            Pac=-1;
            ticker=ticker+1
            stuck=stuck+1;
            idx= ShootPointRND(path); % Choose shooting point
            PsiPert=DisplacementOp((1i*VarMomentum*randn)/sqrt(2),a)*Psi(:,idx);
            PsiPert=PsiPert/norm(PsiPert);
            [newPath, newPsi, idx_new,check,~] =SSEDynamicsTISLeft(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2));
            MP=min(newPath);
            if check==1&&MP<=XValsSSELeft(window)
                Pac=PaccSSETIS(newPsi,Psi,idx,idx_new,ThermalState,HhatGamma,Lhat,hbar,dt);
                if rand <Pac
                    path=newPath;
                    Psi=newPsi;
                    accepted=accepted+1;
                    stuck=0;
                end
            end
            MP=min(path);
            fprintf('SSE Left Window: %5.d, Position %5.4g,  Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,            Pac: %15.2g,   Stuck: %5.d\n', window,XValsSSELeft(window),(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
            if accepted>n_equilib&&mod(ticker,tickrate)==0
                TempMP=[TempMP, MP];
                for n=1:length(MP)
                    bump = abs(customHeavisideLeft(binPositionsLeft - XValsSSELeft(window)) - customHeavisideLeft(binPositionsLeft - MP(n)));
                    TempHist=TempHist+bump;
                end
                TotalTime=TotalTime+dt*length(path);
            end
        end
        MPArrayTemp{copy}=TempMP;
        CrossingHistogramsSSELeft(:,copy,window)=TempHist/length(TempMP);
        TimeTrackerSSELeft(copy,window)=TotalTime/length(TempMP);
    end
    MPArrayTemp=horzcat(MPArrayTemp{:});
    MPArrayTemp=sort(MPArrayTemp(:),"descend");
    MPArraySSELeftWindow{window}=MPArrayTemp;
    XValsSSELeft(window+1)=min(roundToNearest(MPArrayTemp(round((1-Prob)*length(MPArrayTemp))),0.0125),XValsSSELeft(window)-0.0125); % Determines the next window position.
    successfracSSELeft(window)=length(find(MPArrayTemp<=XValsSSELeft(window+1)))/length(MPArrayTemp);
    if window==2
        converged=true;
    else
        CrossingHistogramsSSELeft(:,:,window+1)=zeros(nBinsLeft,nCopy);
        TimeTrackerSSELeft=[TimeTrackerSSELeft, zeros(nCopy,1)];
        fprintf('Window: %5.1g finished,     Next window: %5.2f,    successfrac: %15.3g\n', window,XValsSSELeft(window+1), successfracSSELeft(window))
        window=window+1;
        save(filenameWorkSpace,'-v7.3');
    end
end
fprintf('Window: %5.1g finished, successfrac: %15.3g\n', window, successfracSSELeft(window))
clear MPArrayTemp TempMP TempHist
CrossingHistogramsSSELeft=squeeze(mean(CrossingHistogramsSSELeft,2)).';
[adjusted_weightsSSELeft,unbiasedhistogramSSELeft] =adjust_histogram_weightsTIS(CrossingHistogramsSSELeft,abs(XValsSSELeft),abs(binPositionsLeft),successfracSSELeft);
AvgTMinusSSE=sum(repmat(adjusted_weightsSSELeft.',nCopy,1).*TimeTrackerSSELeft,2);
% AvgTMinusSSE=TimeTrackerSSELeft(:,1);

toc

save(filenameWorkSpace,'-v7.3');
%% Sampling  SSE Right TIS Variable Windows
tic
clear XValsSSE path Psi
XValsSSE=A(2);
window=1;
binPositions = A(2):0.025:B(1);
nBins = length(binPositions);
CrossingHistogramsSSE=zeros(nBins,nCopy,1);
successfracSSE=0;
TimeTrackerSSE=zeros(nCopy,1);
converged=false;
TransitionTimeSSE=zeros(nCopy,1);

%% Loop
while converged==false
    MPArrayTemp=cell(1,nCopy);
    parfor copy = 1:nCopy
        rng(100*copy + sum(clock)*1000);  % Set unique seed based on clock and loop index
        accepted = 0;
        success=0;
        TotalTime=0;
        check=0;
        path=nan;
        Psi=nan;
        stuck=0;
        TempHist=zeros(1,nBins);
        TempMP=[];
        ticker=0;
        while check==0
            PsiPert= DisplacementOpG((XValsSSE(window)+0.1+1i*VarMomentum*randn)/sqrt(2),a,PsiGround);
            [path, Psi, ~,check,~,~] =SSEDynamicsTISRight(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2),B(1));
        end
        while (accepted-n_equilib)< nPathsPerInterface
            Pac=-1;
            ticker=ticker+1
            stuck=stuck+1;
            idx= ShootPointRND(path); % Choose shooting point
            PsiPert=DisplacementOp((1i*VarMomentum*randn)/sqrt(2),a)*Psi(:,idx);
            PsiPert=PsiPert/norm(PsiPert);
            [newPath, newPsi, idx_new,check,~,~] =SSEDynamicsTISRight(PsiPert,dt,tf,HhatGamma,Lhat,Xhat,hbar,A(2),B(1));
            MP=max(newPath);
            if check==1&&MP>=XValsSSE(window)
                Pac=PaccSSETIS(newPsi,Psi,idx,idx_new,ThermalState,HhatGamma,Lhat,hbar,dt);
                if rand <Pac
                    path=newPath;
                    Psi=newPsi;
                    accepted=accepted+1;
                    stuck=0;
                end
            end
            MP=max(path);
            fprintf('SSE Right Window: %5.d, Position %5.4g,  Progress: %5.2f, TAvg: %5.2f,   Steps: %15.3g,            Pac: %15.2g,   Stuck: %5.d\n', window,XValsSSE(window),(accepted-n_equilib)/nPathsPerInterface, TotalTime/length(TempMP),length(path), Pac,stuck)
            if accepted>n_equilib&&mod(ticker,tickrate)==0
                TempMP=[TempMP, MP];
                for n=1:length(MP)
                    bump = customHeaviside(binPositions - XValsSSE(window)) - customHeaviside(binPositions - MP(n));
                    TempHist=TempHist+bump;
                end
                TotalTime=TotalTime+dt*length(path);
            end
        end
        MPArrayTemp{copy}=TempMP;
        CrossingHistogramsSSE(:,copy,window)=TempHist/length(TempMP);
        TimeTrackerSSE(copy,window)=TotalTime/length(TempMP);
    end
    MPArrayTemp=horzcat(MPArrayTemp{:});
    MPArrayTemp=sort(MPArrayTemp(:));
    MPArraySSEWindow{window}=MPArrayTemp;
    XValsSSE(window+1)=max(roundToNearest(MPArrayTemp(round((1-Prob)*length(MPArrayTemp))),0.025),XValsSSE(window)+0.025); % Determines the next window position.
    if XValsSSE(window+1)>=0.5
        converged=true;
        XValsSSE(window+1)=B(1);
        successfracSSE(window)=length(find(MPArrayTemp>=XValsSSE(window+1)))/length(MPArrayTemp);
    else
        CrossingHistogramsSSE(:,:,window+1)=zeros(nBins,nCopy);
        TimeTrackerSSE=[TimeTrackerSSE, zeros(nCopy,1)];
        successfracSSE(window)=length(find(MPArrayTemp>=XValsSSE(window+1)))/length(MPArrayTemp);
        fprintf('Window: %5.1g finished,     Next window: %5.2f,    successfrac: %15.3g\n', window,XValsSSE(window+1), successfracSSE(window))
        window=window+1;
        save(filenameWorkSpace,'-v7.3');
    end
end
fprintf('Window: %5.1g finished, successfrac: %15.3g\n', window, successfracSSE(window))
clear MPArrayTemp TempMP TempHist
CrossingHistogramsSSE=squeeze(mean(CrossingHistogramsSSE,2)).';
[adjusted_weightsSSE,unbiasedhistogramSSE] =adjust_histogram_weightsTIS(CrossingHistogramsSSE,XValsSSE,binPositions,successfracSSE);
AvgTPlusSSE=sum(repmat(adjusted_weightsSSE.',nCopy,1).*TimeTrackerSSE,2);
% AvgTPlusSSE=TimeTrackerSSE(:,1);

AvgTLoopHistSSE=AvgTPlusSSE+AvgTMinusSSE;

winLength=length(successfracSSE);
PTranSSE=prod(successfracSSE);
fprintf('PTranSSE %15.3g\n', PTranSSE)
N_TriesSSE=1./prod(successfracSSE,2);
T_rateSSEHistLoop=1/mean((N_TriesSSE).*AvgTLoopHistSSE)
save(filenameWorkSpace,'-v7.3');
toc

%% transition rate
clear PsiInSSE pathInSSE pathInLangevin PsiInLangevin PsiInSSELeft pathInSSELeft pathInLangevinLeft PsiInLangevinLeft Psi path

N_TriesLangevin=1./prod(successfracLangevin,2);
T_rateLangevinHistLoop=1/mean((N_TriesLangevin).*AvgTLoopHistLangevin)
N_TriesSSE=1./prod(successfracSSE,2);
T_rateSSEHistLoop=1/mean((N_TriesSSE).*AvgTLoopHistSSE)

save(filenameWorkSpace,'-v7.3');
%%
set(0, 'DefaultTextInterpreter', 'latex')
CrossingHistogramsLangevin(CrossingHistogramsLangevin == 0) = 1;
CrossingHistogramsSSE(CrossingHistogramsSSE == 0) = 1;

f=figure
hold on
for n=1:length(successfracLangevin)
    plot(binPositions,CrossingHistogramsLangevin(n,:),'--')
end
% plot(binPositions,Prob*ones(nBins),'-r',LineWidth=1)
xlim([A(2) B(1)])
xticks(linspace(A(2), B(1),5))
yticks(linspace(0, 1,6))

ylim([0, 1])
xlabel('$\lambda$')
ylabel('$P(\lambda|\lambda_i)$')
box on
filename = sprintf('CrossingPlotLangevin.pdf');
ax = gca;
ax.XTickLabelRotation = 0; % Add this line
p.LineWidth = 1;
f.Position = [0 0 860 500]; % Set the width and height of the figure
ax.FontSize = 25;  % Set to desired font size
% exportgraphics(f,filename,'ContentType','vector')

f=figure
hold on
for n=1:length(successfracSSE)
    plot(binPositions,CrossingHistogramsSSE(n,:),'--')
end
% plot(binPositions,Prob*ones(nBins),'-r',LineWidth=1)
xlim([A(2) B(1)])
xticks(linspace(A(2), B(1),5))
yticks(linspace(0, 1,6))
ylim([0, 1])
xlabel('$\lambda$')
ylabel('$P(\lambda|\lambda_i)$')
box on
filename = sprintf('CrossingPlotSSE.pdf');
ax = gca;
ax.XTickLabelRotation = 0; % Add this line
p.LineWidth = 1;
f.Position = [0 0 860 500]; % Set the width and height of the figure
ax.FontSize = 25;  % Set to desired font size
% exportgraphics(f,filename,'ContentType','vector')

%% functions

function result = DisplacementOp(alphaVec,a)
N = length(alphaVec);
bSize=size(a,1);
result = zeros(bSize, bSize, N);
for n = 1:N
    alpha = alphaVec(n);
    result(:,:,n) = expm(alpha*a' - conj(alpha)*a);
end
end

function result = DisplacementOpG(alphaVec,a,PsiGround)
N = length(alphaVec);
bSize=size(a,1);
result = zeros(bSize, N);
for n = 1:N
    alpha = alphaVec(n);
    result(:,n) = expm(alpha*a' - conj(alpha)*a)*PsiGround;
end
end


function y = customHeaviside(x)
y = double(x >= 0);  % Returns 1 for x >= 0 and 0 otherwise
end

function y = customHeavisideLeft(x)
y = double(x > 0);  % Returns 1 for x >= 0 and 0 otherwise
end

function roundedValue = roundToNearest(value,num)
roundedValue = round(value / num) * num;
end