function U= BM3D_Hybrid(U1,I2,sigma,opts)

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


parfor ii = 1:numel(opts.matchSpins)
    % second estimate: re-match blocks and empirical Wiener filter in
    % transform domain (2D DCT + Haar)
    [S,V1,V2] = matchBlocksWie(...
        circshift(U1,[opts.matchSpins(ii),opts.matchSpins(ii)]),...
        circshift(I2,[opts.matchSpins(ii),opts.matchSpins(ii)]),...
        sigma,opts);
    [U(:,:,ii),out] = denoiseAndAggregateWie(V1,V2,S,sigma,opts);
    U(:,:,ii) = circshift(U(:,:,ii),[-opts.matchSpins(ii),-opts.matchSpins(ii)]);
end
% trim images back to original dimensions and average over the spins
% out.U1 = U1(1:m,1:n);
U = mean(U,3); % final estimate
