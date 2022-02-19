function [U] = LTBM3Dparallel(I,sigma,opts)

if ~isfield(opts,'blockSize'), opts.blockSize = 64; end
if ~isfield(opts,'numMax'), opts.numMax = 64; end
if ~isfield(opts,'numMin'), opts.numMin = 32; end
if ~isfield(opts,'wname'), opts.wname = 'sym'; end
if ~isfield(opts,'wnamez'), opts.wnamez = 'db'; end
if ~isfield(opts,'order'), opts.order = 4; end
if ~isfield(opts,'orderz'), opts.orderz = 1; end
if ~isfield(opts,'levels'), opts.levels = 3; end
if ~isfield(opts,'cycleSpin'), opts.cycleSpin = 1; end
if ~isfield(opts,'wiener'), opts.wiener = true; end
if ~isfield(opts,'matchSpins'), opts.matchSpins = 0; end
if ~isfield(opts,'tauMode'), opts.tauMode = 2; end

% reflect the image so that its dimension is a power of two
% it likely only needs to be divisible by the blocksize, but this is
% suitable for now
[m,n] = size(I);
M = 2^ceil(log(m)/log(2));
N = 2^ceil(log(n)/log(2));
I2 = zeros(M,N);
I2(1:m,1:n) = I;
I2(m+1:M,1:n) = flipud(I(m-(M-m)+1:m,1:n));
I2(1:m,n+1:N) = fliplr(I(1:m,n-(N-n)+1:n));
I2(m+1:M,n+1:N) = fliplr(flipud(I(m-(M-m)+1:m,n-(N-n)+1:n)));
U = zeros(M,N,numel(opts.matchSpins));

opts.filtType = 'ht';
% loop of "match spins" where the algorithm is re-evaluated after circle
% shifting the image, so that brand new set of reference blocks are used
parfor ii = 1:numel(opts.matchSpins)
    % match blocks and denoise
      %   I2 = circshift(I2,[opts.matchSpins(ii),opts.matchSpins(ii)]);
        S = matchBlocksGPU(I2,sigma,opts);
        [U(:,:,ii),~,~] = denoiseAndAggregateGPU(...
            circshift(I2,[opts.matchSpins(ii),opts.matchSpins(ii)]),S,sigma,opts);
        U(:,:,ii) = circshift(U(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
  %      I2 = circshift(I2,[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
end
% trim the image back to original dimensions and average over the spins
U = mean(U(1:m,1:n,:),3);
