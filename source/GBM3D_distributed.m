function [U,out] = GBM3D_distributed(I,sigma,opts)

% G-BM3D: a fast variation of BM3D with fast block matching, global
% wavelet thresholding, and GPU capabilities.
% This version is parallelized further by distributing the image into
% patches and processing in parallel

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
[m,n] = size(I);

% if image is too small, run standard algorithm
if m<256 + opts.blockSize || n<256 + opts.blockSize 
    [U,out] = GBM3D(I,sigma,opts);
    return;
end

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

% set up the image distribution for parallel processing
[R,d0,K] = getBM3D_distribution(I,opts.matchSpins,opts.blockSize);
[RW,d0,K] = getBM3D_distribution(I,opts.matchSpinsWie,opts.blockSize);

% set up cell variables to store processed image subsets
U1 = cell(size(R,1),1);
W = U1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIRST ESTIMATE: wavelet hard thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over each case (image patch + translation)
parfor ii = 1:size(R,1)
    indY = R(ii,1):R(ii,2);
    indX = R(ii,3):R(ii,4);
    Ilocal = circshift(I(indY,indX),[R(ii,5),R(ii,5)]);
    [S,Volume] = matchBlocksCPU(Ilocal,sigma,opts);
    [U1{ii},W{ii}] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts);
    U1{ii} = circshift(U1{ii},[-R(ii,5),-R(ii,5)]);
    W{ii} = circshift(W{ii},[-R(ii,5),-R(ii,5)]);
end


% stitch all of the cases back together into U2 to obtain first estimate
U2 = zeros(m,n);
mWindow = makeInterpWindow(d0,d0,K); % window used for most cases
normMatrix = zeros(m,n); % normalization matrix with weights
for ii = 1:size(R,1)
    indY = R(ii,1):R(ii,2);
    indX = R(ii,3):R(ii,4);
    if numel(indX)~=d0 || numel(indY)~=d0
        tWindow = makeInterpWindow(numel(indY),numel(indX),K);
    else
        tWindow = mWindow;
    end
    U2(indY,indX) = U2(indY,indX) + U1{ii}.*tWindow;
    normMatrix(indY,indX) = normMatrix(indY,indX) + tWindow.*W{ii};
end
U2 = U2./normMatrix; % first estimate
out.U2 = U2;
if ~opts.Wiener
    U = U2(1:m,1:n);
    return;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      SECOND ESTIMATE: empirical Wiener filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over each case (image patch + translation)
parfor ii = 1:size(RW,1)
    indY = RW(ii,1):RW(ii,2);
    indX = RW(ii,3):RW(ii,4);
    Ilocal1 = circshift(I(indY,indX),[RW(ii,5),RW(ii,5)]);
    Ilocal2 = circshift(U2(indY,indX),[RW(ii,5),RW(ii,5)]);

    [S,V1,V2] = matchBlocksWie(Ilocal2,Ilocal1,sigma,opts);
    [U1{ii},W{ii}] = denoiseAndAggregateWieFast(V1,V2,S,sigmaWie,opts);
    U1{ii} = circshift(U1{ii},[-RW(ii,5),-RW(ii,5)]);
    W{ii} = circshift(W{ii},[-RW(ii,5),-RW(ii,5)]);
end

% stitch all of the cases back together to obtain final estimate
U2 = zeros(m,n);
normMatrix = zeros(m,n); % normalization matrix with weights
for ii = 1:size(RW,1)
    indY = RW(ii,1):RW(ii,2);
    indX = RW(ii,3):RW(ii,4);
    if numel(indX)~=d0 || numel(indY)~=d0
        tWindow = makeInterpWindow(numel(indY),numel(indX),K);
    else
        tWindow = mWindow;
    end
    U2(indY,indX) = U2(indY,indX) + U1{ii}.*tWindow;
    normMatrix(indY,indX) = normMatrix(indY,indX) + tWindow.*W{ii};
end
U = U2./normMatrix; % final estimate


    