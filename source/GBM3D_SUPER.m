function [U,out] = GBM3D_SUPER(I,K,hhat,sigma,opts)

% G-BM3D: a fast variation of BM3D with fast block matching, global
% wavelet thresholding, and GPU capabilities

% Reference: Toby Sanders and Sean Larkin. "New computational techniques
% for a faster variation of BM3D." Submitted in 2021.


% Written by Toby Sanders @Lickenbrock Tech
% Date: 5-9-21


% check for the selected profile and set all GBM3D options
if nargin<3, opts.profile = 'accuracy'; end
if isfield(opts,'profile')
    opts = setBM3Dopts(opts.profile); 
end
opts = checkBM3Dopts(opts);


[m,n] = size(I);
[m2,n2] = size(hhat);
if m2~=K*m || n2~=K*n
    error('PSF dimensions should match upsampling resolution');
end
if abs(hhat(round(m2/2),round(n2/2)))>1e-2
    useSolver = true;
else
    useSolver = false;
end
useSolver = true;


Istretch = zeros(m*K,n*K);
Istretch(1:K:end,1:K:end) = I;
FI = fft2(Istretch);
hhat2 = abs(hhat).^2;
M = ceil(m*K/opts.blockSize)*opts.blockSize;
N = ceil(n*K/opts.blockSize)*opts.blockSize;
V = my_Fourier_filters(1,1,m*K,n*K,1);
if numel(sigma)==1
    sigma = ones(m*K,n*K)*sigma^2; 
end
alphaRI = 250*sigma;
alphaRWI = alphaRI*10*4;

g = zeros(m*K,n*K);
g([1:K],[1:K]) = 1/K^2;
g = fraccircshift(g,[-K/2 + 1/2, -K/2 + 1/2]);
ghat = fft2(g);
ghat2 = abs(ghat).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FIRST WIENER DECON. FILTER ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filtRI = conj(hhat).*conj(ghat)./(hhat2.*ghat2/K^2 + alphaRI.*V);
out.recWie = real(ifft2(filtRI.*FI));
if useSolver
    fudgeFactor = 1e-1;
    regMatrix = fudgeFactor*V.*alphaRI;
    [out.recWie,out.tik1] = Tikhonov_SUPER(I,K,hhat,regMatrix);
end
sigmaPSD = sigma.*abs(filtRI).^2;
tau = getColoredTau(opts.levels,sigmaPSD);

% get 3D wavelet filters (no longer any inputs, default is always
% biorthogonal wavelets order 1.5 crossed with Haar wavelet)
[Psi,Psi2] = getWaveFilters3Dbior;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIRST Denoise ESTIMATE: wavelet hard thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop of "match spins" where the algorithm is re-evaluated after circle
% shifting the image, so that brand new set of reference blocks are used 
matchSpins = opts.matchSpins;
Itmp = myMirrorPadImage(out.recWie,M,N); % pad image so its dims are divisible by blocksize
U = zeros(M,N,numel(matchSpins));
W = U;U1 = U;
parfor ii = 1:numel(matchSpins)
    [S,Volume] = matchBlocksCPU(circshift(Itmp,[matchSpins(ii),matchSpins(ii)]),sigma,opts);
    [U1(:,:,ii),W(:,:,ii)] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts);
    U1(:,:,ii) = circshift(U1(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
    W(:,:,ii) = circshift(W(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
end
% average over the translations
U1 = sum(U1,3)./sum(W,3); % first estimate

if ~opts.Wiener % if not performing second step
    out.U1 = U1(1:m*K,1:n*K);
    U = U1(1:m*K,1:n*K); % final estimate
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SECOND WIENER DECON. FILTER ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wiener_Pilot = abs(fft2(U1(1:m*K,1:n*K))).^2;
filtRWI = conj(hhat).*conj(ghat).*Wiener_Pilot./(Wiener_Pilot.*hhat2.*ghat2/K^2 + alphaRWI);
if useSolver
    fudgeFactor = 1e-4;
    regMatrix = fudgeFactor*alphaRWI./Wiener_Pilot;
    [out.recWie2,out.tik2] = Tikhonov_SUPER(I,K,hhat,regMatrix);
end
% out.recWie2 = real(ifft2(FI.*filtRWI));
sigmaPSD = sigma.*abs(filtRWI).^2;
sigmaWie = getColoredWie(sigmaPSD,opts.blockSizeWie);


Itmp = myMirrorPadImage(out.recWie2,M,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      SECOND DENOISE ESTIMATE: empirical Wiener filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor ii = 1:numel(matchSpins)
    [S,V1,V2] = matchBlocksWie(...
        circshift(U1,[matchSpins(ii),matchSpins(ii)]),...
        circshift(Itmp,[matchSpins(ii),matchSpins(ii)]),...
        sigma,opts);
    [U(:,:,ii),W(:,:,ii)] = denoiseAndAggregateWieFast(V1,V2,S,sigmaWie,opts);
    U(:,:,ii) = circshift(U(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
    W(:,:,ii) = circshift(W(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
end
% trim images back to original dimensions and average over the spins
out.U1 = U1(1:m*K,1:n*K);
U = sum(U,3)./sum(W,3);
U = U(1:m*K,1:n*K);  % final estimate

