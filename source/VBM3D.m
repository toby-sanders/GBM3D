function [Y,out] = VBM3D(I,sigma,opts)

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
if ~isfield(opts,'wiener'), opts.wiener = false; end
if ~isfield(opts,'matchSpins'), opts.matchSpins = [0,opts.blockSize/2]; end
if ~isfield(opts,'tauMode'), opts.tauMode = 2; end
if ~isfield(opts,'nF'), opts.nF = 11; end

% reflect the image so that its dimension is a power of two
% it likely only needs to be divisible by the blocksize, but this is
% suitable for now
[m,n,K] = size(I);
% M = 2^ceil(log(m)/log(2));
% N = 2^ceil(log(n)/log(2));
M = ceil(m/opts.blockSize)*opts.blockSize;
N = ceil(n/opts.blockSize)*opts.blockSize;
I2 = zeros(M,N,K);
I2(1:m,1:n,:) = I;
I2(m+1:M,1:n,:) = flipud(I(m-(M-m)+1:m,1:n,:));
I2(1:m,n+1:N,:) = fliplr(I(1:m,n-(N-n)+1:n,:));
I2(m+1:M,n+1:N,:) = fliplr(flipud(I(m-(M-m)+1:m,n-(N-n)+1:n,:)));
U = zeros(M,N,opts.nF,numel(opts.matchSpins));


% loop of "match spins" where the algorithm is re-evaluated after circle
% shifting the image, so that brand new set of reference blocks are used

% aggregation will later be stored in these variables
Y = zeros(M,N,K,numel(opts.matchSpins),'single');
nn= Y(:,:,:,1);
f1 = (opts.nF+1)/2;
f2 = K-(opts.nF-1)/2;
for ii = 1:numel(opts.matchSpins)
    nn(:) = 0;
    for jj = f1:f2
        ind = jj-(opts.nF-1)/2:jj+(opts.nF-1)/2;
        I3 = I2(:,:,ind);
    % match blocks and denoise
        [S,Volume] = matchBlocksV(circshift(I3,[opts.matchSpins(ii),opts.matchSpins(ii),0]),sigma,opts);
        [U,out] = denoiseAndAggregateV(Y,Volume,...
            circshift(I3,[opts.matchSpins(ii),opts.matchSpins(ii),0]),S,sigma,opts);     
        Y(:,:,ind,ii) = Y(:,:,ind,ii) + U;
        nn(:,:,ind) = nn(:,:,ind) + out.K;
    end
    Y(:,:,:,ii) = Y(:,:,:,ii)./nn;
    % U(:,:,:,ii) = U(:,:,:,ii)./nn;
    Y(:,:,:,ii) = circshift(Y(:,:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii),0]);    
end
% trim the image back to original dimensions and average over the spins
Y = mean(Y(1:m,1:n,f1:f2,:),4);
