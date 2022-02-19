function [U,out] = LTBM3D(I,sigma,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

% a GPU variation of BM3D 

if ~isfield(opts,'blockSize'), opts.blockSize = 32; end
if ~isfield(opts,'numMax'), opts.numMax = 16; end
if ~isfield(opts,'numMin'), opts.numMin = 16; end
if ~isfield(opts,'wname'), opts.wname = 'sym'; end
if ~isfield(opts,'wnamez'), opts.wnamez = 'db'; end
if ~isfield(opts,'order'), opts.order = 4; end
if ~isfield(opts,'orderz'), opts.orderz = 1; end
if ~isfield(opts,'levels'), opts.levels = 3; end
if ~isfield(opts,'cycleSpin'), opts.cycleSpin = 2; end
if ~isfield(opts,'matchSpins'), opts.matchSpins = [0,opts.blockSize/2]; end
if ~isfield(opts,'tauMode'), opts.tauMode = 2; end
if ~isfield(opts,'filtType'), opts.filtType = 'ht'; end
if ~isfield(opts,'blockSizeWie'), opts.blockSizeWie = 16; end
if ~isfield(opts,'matchSize'), opts.matchSize = opts.blockSize; end
if ~isfield(opts,'Wiener'), opts.Wiener = true; end



% pad image so its dimensions are divisible by the blocksize
[m,n] = size(I);
M = ceil(m/opts.blockSize)*opts.blockSize;
N = ceil(n/opts.blockSize)*opts.blockSize;

% reflect the image to pad to new dimension
I2 = zeros(M,N);
I2(1:m,1:n) = I;
I2(m+1:M,1:n) = flipud(I(m-(M-m)+1:m,1:n));
I2(1:m,n+1:N) = fliplr(I(1:m,n-(N-n)+1:n));
I2(m+1:M,n+1:N) = fliplr(flipud(I(m-(M-m)+1:m,n-(N-n)+1:n)));
U = zeros(M,N,numel(opts.matchSpins));
W = U;
W1 = U;
U1 = U;

% precompute threshold levels tau and wavelet filters
if size(sigma,1)==1 
    tau = sigma*getTau(opts.levels,opts.tauMode);
else
    tau = getColoredTau(opts.levels,opts.order,sigma);
    opts.Wiener = false;
end
if strcmp(opts.wname,'bior')
    [Psi,Psi2] = getWaveFilters3Dbior(opts.wnamez,opts.order,opts.orderz);
else
    Psi = getWaveFilters3D(opts.wname,opts.wnamez,opts.order,opts.orderz);
    Psi2 = Psi;
end

% loop of "match spins" where the algorithm is re-evaluated after circle
% shifting the image, so that brand new set of reference blocks are used 
for ii = 1:numel(opts.matchSpins)
    % first etimate: match blocks and wavelet hard thresholding
    [S,Volume] = matchBlocksCPU(circshift(I2,[opts.matchSpins(ii),opts.matchSpins(ii)]),sigma,opts);
    [U1(:,:,ii),W1(:,:,ii)] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts);
    U1(:,:,ii) = circshift(U1(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
    W1(:,:,ii) = circshift(W1(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
end
% average over the spins
U1 = mean(U1./W1,3); % first estimate

% second estimate: re-match blocks and empirical Wiener filter in
% transform domain (2D DCT + Haar)
if opts.Wiener
    parfor ii = 1:numel(opts.matchSpins)
        [S,V1,V2] = matchBlocksWie(...
            circshift(U1,[opts.matchSpins(ii),opts.matchSpins(ii)]),...
            circshift(I2,[opts.matchSpins(ii),opts.matchSpins(ii)]),...
            sigma,opts);
        [U(:,:,ii),W(:,:,ii)] = denoiseAndAggregateWie(V1,V2,S,sigma,opts);
        U(:,:,ii) = circshift(U(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
        W(:,:,ii) = circshift(W(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
    end
    % trim images back to original dimensions and average over the spins
    out.U1 = U1(1:m,1:n);
    U = mean(U./W,3);
    U = U(1:m,1:n);
    % U = mean(U(1:m,1:n,:),3); % final estimate
else
    out.U1 = U1(1:m,1:n);
    U = U1(1:m,1:n); % final estimate
end
