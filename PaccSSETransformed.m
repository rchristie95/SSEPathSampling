function Pacc=PaccSSETransformed(SPsiN,PsiO,idx,ThermalState,Hhat,Lhat,hbar,dt,choice)
bSize=size(PsiO,1);
N=size(PsiO,2);
Parity=((-ones(1,bSize)).^(0:bSize-1)).';
Parity=repmat(Parity,1,N); % pagevector parity operator fock basis

if choice==2 || choice==4
    idxPrime=N-(idx-1);
else
    idxPrime=idx;
end

ChoiceString={'S','T','P','PT'};
functionName = strcat(ChoiceString{choice}, '_Reversal');
PsiN=feval(functionName, SPsiN,Parity);
PsiSO=feval(functionName, PsiO,Parity);

P0=(SPsiN(:,1)'*ThermalState*SPsiN(:,1))/(PsiO(:,1)'*ThermalState*PsiO(:,1));
PsiNB=T_Reversal(PsiN(:,1:idx),Parity);
PsiSOB=T_Reversal(PsiSO(:,1:idxPrime),Parity);
S_SN=zeros(1,N);
S_SO=zeros(1,N);
S_N=zeros(1,N);
S_O=zeros(1,N);

for n=1:N-1
    Psi0=SPsiN(:,n)/norm(SPsiN(:,n));
    lev=SPsiN(:,n)'*Lhat*SPsiN(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((SPsiN(:,n+1)-Psi0)/dt-mu);
    S_SN(n)=(top'*top)/2/(sigma'*sigma)^2;

    Psi0=PsiO(:,n)/norm(PsiO(:,n));
    lev=PsiO(:,n)'*Lhat*PsiO(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((PsiO(:,n+1)-Psi0)/dt-mu);
    S_O(n)=(top'*top)/2/(sigma'*sigma)^2;
end

for n=1:idx-1
    Psi0=PsiNB(:,n)/norm(PsiNB(:,n));
    lev=PsiNB(:,n)'*Lhat*PsiNB(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((PsiNB(:,n+1)-Psi0)/dt-mu);
    S_N(n)=(top'*top)/2/(sigma'*sigma)^2;
end

for n=idx:N-1
    Psi0=PsiN(:,n)/norm(PsiN(:,n));
    lev=PsiN(:,n)'*Lhat*PsiN(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((PsiN(:,n+1)-Psi0)/dt-mu);
    S_N(n)=(top'*top)/2/(sigma'*sigma)^2;
end

for n=1:idxPrime-1
    Psi0=PsiSOB(:,n)/norm(PsiSOB(:,n));
    lev=PsiSOB(:,n)'*Lhat*PsiSOB(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((PsiSOB(:,n+1)-Psi0)/dt-mu);
    S_SO(n)=(top'*top)/2/(sigma'*sigma)^2;
end

for n=idxPrime:N-1
    Psi0=PsiSO(:,n)/norm(PsiSO(:,n));
    lev=PsiSO(:,n)'*Lhat*PsiSO(:,n)*eye(bSize);
    sigma=(1/sqrt(hbar))*(Lhat-lev)*Psi0;
    mu=hbar^(-1)*(-1i*Hhat-0.5*(Lhat'*Lhat)+lev'*Lhat-0.5*(lev'*lev))*Psi0;
    top=(sigma)'*((PsiSO(:,n+1)-Psi0)/dt-mu);
    S_SO(n)=(top'*top)/2/(sigma'*sigma)^2;
end

S=dt*sum(S_SN-S_O+S_SO-S_N);
Pacc=real(P0*exp(-S));
end
%% functions

function [outPsi]= S_Reversal(Psi,Parity)
outPsi=Psi;
end

function [outPsi]= T_Reversal(Psi,Parity)
outPsi=conj(flip(Psi,2));
end

function [outPsi]= P_Reversal(Psi,Parity)
outPsi=Parity.*Psi;
end

function [outPsi]= PT_Reversal(Psi,Parity)
outPsi=conj(flip(Parity.*Psi,2));
end