function out = FockWigner(PsiIn,x,p,hbar,HVector)
%MYWIGNER: Calculates the Wigner distribution from a column vector
%
%	W  = mywigner(Ex)
%
%	W  = output Wigner distribution
%	Ex = Input electric field (MUST be a column vector
%
%	Notes:
%		W = Int(-inf..inf){E(x+y)E(x-y)exp[2ixy]}
%
%		E(x+y) & E(x-y) are calculated via a FFT (fast Fourier transform) using the
%		shift theorem. The integration is performed via a FFT. Thus it is important
%		for the data to satisfy the sampling theorem:
%		dy = 2*pi/X			X = span of all x-values	dy = y resolution
%		dx = 2*pi/Y			Y = span of all y-values	dx = x resolution
%		The data must be completely contained within the range x(0)..x(N-1) &
%		y(0)..y(N-1) (i.e. the function must fall to zero within this range).
%
%	v1.0
%
%	Currently waiting for update:
%		Remove the fft/ifft by performing this inside the last function calls
%		Allow an arbitrary output resolution
%		Allow an input vector for x (and possibly y).




N=size(PsiIn,2);
Nx = length(x); %   Get length of vector
x2=x.';
p2=p.';
x2=ifftshift(x2);
% W=zeros(Nx,Nx,N);
W=zeros(Nx,Nx);
% Wig=zeros(Nx,Nx,N);
for n=1:N
    Ex=HVector*PsiIn(:,n);
    EX1 = ifft( (fft(Ex)*ones(1,Nx)).*exp( 1i*x2*p2.'/2/hbar ));					%   +ve shift
    EX2 = ifft( (fft(Ex)*ones(1,Nx)).*exp( -1i*x2*p2.'/2/hbar ));					%   -ve shift
    tempWig=(1/2/pi/hbar)*real(fftshift(fft(fftshift(EX1.*conj(EX2), 2), [], 2), 2))';
    % Wig(:,:,n) =tempWig;		%   Wigner function
    W= W+tempWig;		%   Wigner function
end

% out=W;
out=W/N;


%% Wigner Videos
% figure
% Xview=6;
% [x1, p1]=meshgrid(x,p);
% 
% f=surf(x1,p1,Wig(:,:,1));
% shading interp
% % title('Psi1')
% xlim([-Xview Xview])
% ylim([-Xview Xview])
% 
% 
% xlabel('$x$',Interpreter='latex')
% ylabel('$p$',Interpreter='latex')
% view(2)
% caxis([-0.05 0.3])
% colormap(jet)
% xticks(-Xview:Xview)
% 
% for n=2:length(Wig(1,1,:))
%     set(f,'ZData',Wig(:,:,n));
%     drawnow
% end

end
