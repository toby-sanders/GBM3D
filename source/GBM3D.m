function [U,out] = GBM3D(I,sigma,opts)

% G-BM3D: a fast variation of BM3D with fast block matching, global
% wavelet thresholding, and GPU capabilities

% Reference: Toby Sanders and Sean Larkin. "New computational techniques
% for a faster variation of BM3D." Submitted in 2021.


% Written by Toby Sanders @Lickenbrock Tech
% 3-17-2021


% check for the selected profile and set all GBM3D options
if nargin<3, opts.profile = 'default'; end
if isfield(opts,'profile')
    opts = setBM3Dopts(opts.profile); 
end
opts = checkBM3Dopts(opts);

% pad image so its dimensions are divisible by the blocksize
[m,n] = size(I);
M = ceil(m/opts.blockSize)*opts.blockSize;
N = ceil(n/opts.blockSize)*opts.blockSize;
I = myMirrorPadImage(I,M,N);
U = zeros(M,N,numel(opts.matchSpins));
W = U;U1 = U;

% precompute threshold levels tau
if size(sigma,1)==1 
    tau = sigma*getTau(opts.levels,opts.tauMode);
    sigmaWie = sigma;
else
    tau = getColoredTau(opts.levels,sigma);
    sigmaWie = getColoredWie(sigma,opts.blockSizeWie);
end

% get 3D wavelet filters (no longer any inputs, default is always
% biorthogonal wavelets order 1.5 crossed with Haar wavelet)
[Psi,Psi2] = getWaveFilters3Dbior;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIRST ESTIMATE: wavelet hard thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop of "match spins" where the algorithm is re-evaluated after circle
% shifting the image, so that brand new set of reference blocks are used 
matchSpins = opts.matchSpins;
for ii = 1:numel(matchSpins)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      SECOND ESTIMATE: empirical Wiener filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matchSpins = opts.matchSpinsWie;
for ii = 1:numel(matchSpins)
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

