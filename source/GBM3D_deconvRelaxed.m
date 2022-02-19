function [U,out] = GBM3D_deconv(I,hhat,sigma,opts)

% G-BM3D: a fast variation of BM3D with fast block matching, global
% wavelet thresholding, and GPU capabilities

% Reference: Toby Sanders and Sean Larkin. "New computational techniques
% for a faster variation of BM3D." Submitted in 2021.

% Written by Toby Sanders @Lickenbrock Tech
% Date: 5-9-21

% check if BM3D options were passed in
if nargin<3, opts.profile = 'accuracy'; end

% check if alpha regularization parameter is passed in as an option
% if not, use the default tuned alpha based on sigma
% do this before calling "setBM3Dopts" since it will erase the input alpha
[m,n] = size(I);
if numel(sigma)==1
    sigma = ones(m,n)*sigma^2; 
end
if ~isfield(opts,'alpha')
    alphaRI = 27*sigma;
    alphaRWI = alphaRI*10*4;
else
    alphaRI = opts.alpha*ones(m,n);
    alphaRWI = alphaRI*10*4;
end

% check for the selected profile and set all GBM3D options
if isfield(opts,'profile')
    opts = setBM3Dopts(opts.profile); 
end
opts = checkBM3Dopts(opts);


% initialize several variables
FI = fft2(I);
hhat2 = abs(hhat).^2;
M = ceil(m/opts.blockSize)*opts.blockSize;
N = ceil(n/opts.blockSize)*opts.blockSize;
V = my_Fourier_filters(1,1,m,n,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FIRST WIENER DECON. FILTER ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parm.theta = sigma^2*50;parm.order = 1;
% [~,out.ME] = HOTVL2_deblurMF(ifft2(hhat),I,parm);
% alphaRI = out.ME.thetas(end)/2;
filtRI = conj(hhat)./(hhat2 + alphaRI.*V);
out.recWie = real(ifft2(filtRI.*FI));
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
I = myMirrorPadImage(out.recWie,M,N); % pad image so its dims are divisible by blocksize
U = zeros(M,N,numel(matchSpins));
W = U;U1 = U;
parfor ii = 1:numel(matchSpins)
    [S,Volume] = matchBlocksCPU(circshift(I,[matchSpins(ii),matchSpins(ii)]),sigma,opts);
    [U1(:,:,ii),W(:,:,ii)] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts);
    U1(:,:,ii) = circshift(U1(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
    W(:,:,ii) = circshift(W(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
end
% average over the translations
U1 = sum(U1,3)./sum(W,3); % first estimate

if ~opts.Wiener % if not performing second step
    out.U1 = U1(1:m,1:n);
    U = U1(1:m,1:n); % final estimate
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SECOND WIENER DECON. FILTER ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wiener_Pilot = 1./(abs(fft2(U1(1:m,1:n))).^2 + 1e-4);


% relax the regularization by averaging with first step regularizer
gamma = 0.5;
Reg = (gamma*V.*alphaRI + (1-gamma)*Wiener_Pilot.*alphaRWI);
% gamma = 1e-2;
% scl = sum(V(:).*Wiener_Pilot(:))./sum(V(:).^2)/2;
% Reg = (scl*gamma*V + (1-gamma)*Wiener_Pilot).*alphaRWI;


filtRWI = conj(hhat)./(hhat2 + Reg);
out.recWie2 = real(ifft2(FI.*filtRWI));
sigmaPSD = sigma.*abs(filtRWI).^2;
sigmaWie = getColoredWie(sigmaPSD,opts.blockSizeWie);


I = myMirrorPadImage(out.recWie2,M,N);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      SECOND DENOISE ESTIMATE: empirical Wiener filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor ii = 1:numel(matchSpins)
    [S,V1,V2] = matchBlocksWie(...
        circshift(U1,[matchSpins(ii),matchSpins(ii)]),...
        circshift(I,[matchSpins(ii),matchSpins(ii)]),...
        sigma,opts);
    [U(:,:,ii),W(:,:,ii)] = denoiseAndAggregateWieFast(V1,V2,S,sigmaWie,opts);
    U(:,:,ii) = circshift(U(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
    W(:,:,ii) = circshift(W(:,:,ii),[-matchSpins(ii),-matchSpins(ii)]);
end
% trim images back to original dimensions and average over the spins
out.U1 = U1(1:m,1:n);
U = sum(U,3)./sum(W,3);
U = U(1:m,1:n);  % final estimate

