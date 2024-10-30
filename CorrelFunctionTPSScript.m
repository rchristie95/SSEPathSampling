
%% Data Entry
% license('test', 'Distrib_Computing_Toolbox')
% % Manually specify the number of cores requested
% numCores = 16; % Update this with the number of cores requested in your PBS script
% fprintf('Number of allocated cores: %d\n', numCores);
% 
% poolobj = gcp('nocreate');
% if ~isempty(poolobj)
%     delete(poolobj);
% end
% parpool('local', numCores);
% 
% % Display the number of workers in the parallel pool
% currentPool = gcp;
% fprintf('Number of workers in the current parallel pool: %d\n', currentPool.NumWorkers);

%% Main

% % %Quartic
% Xview=4.5;
% h2=0.35; %Quartic
% h4=(h2+0.5)/(Xview^2);
% A=[-3.5;-1];
% B=flip(abs(A));

h2=0.35; %Quartic
h4=0.01;
A=[-6;-4];
B=flip(abs(A));

% nPathsPerWindow=1;
% n_equilib=1;
nPathsPerWindow=2000;
n_equilib=100;
% windowMat=[-6:1.5:4.5;-4:1.5:6.5];
windowMat=[-6.1:1.5:4.4;-4.4:1.5:6.1];

nWindows=size(windowMat,2);
LambdaVals=0.5*(windowMat(1,:)+windowMat(2,:));
% windowMat(1,:)=LambdaVals-WindowWidth/2;
% windowMat(2,:)=LambdaVals+WindowWidth/2;
%% Langevin

LangevinTPSSymmetry(h2,h4,0.25,1.5312/2,LambdaVals,windowMat,A,n_equilib,nPathsPerWindow);
SSETPSSymmetry(h2,h4,0.25,1.5312/2,LambdaVals,windowMat,A,n_equilib,nPathsPerWindow)

% delete(gcp);