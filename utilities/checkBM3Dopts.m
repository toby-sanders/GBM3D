function opts = checkBM3Dopts(opts)

% check that all the options for the GBM3D algorithm have been specified.
% If not set the default recommended values given below. The description of
% each option is given at the bottom

% Written by Toby Sanders @Lickenbrock Tech
% 3-17-2021

if ~isfield(opts,'numMax'), opts.numMax = 16; end
if ~isfield(opts,'numMin'), opts.numMin = 16; end
if ~isfield(opts,'wname'), opts.wname = 'bior'; end
if ~isfield(opts,'wnamez'), opts.wnamez = 'db'; end
if ~isfield(opts,'order'), opts.order = 15; end
if ~isfield(opts,'orderz'), opts.orderz = 1; end
if ~isfield(opts,'levels'), opts.levels = 3; end
if ~isfield(opts,'tauMode'), opts.tauMode = 1; end
if ~isfield(opts,'filtType'), opts.filtType = 'ht'; end
if ~isfield(opts,'blockSize'), opts.blockSize = 16; end
if ~isfield(opts,'blockSizeWie'), opts.blockSizeWie = 8; end
if ~isfield(opts,'matchSize'), opts.matchSize = max(32,opts.blockSize); end
if ~isfield(opts,'cycleSpin'), opts.cycleSpin = 2; end
if ~isfield(opts,'matchSpins'), opts.matchSpins = [0,min(opts.blockSize,opts.blockSizeWie)/2]; end
if ~isfield(opts,'Wiener'), opts.Wiener = true; end

% DESCRIPTION OF EACH OPTION
% -------------------------------------------------------------------
% numMax - maximum number of matched blocks for each reference block
% numMin - minimum number of matched blocks for each reference block
% wname - name of wavelet used for 2D transform
% wnamez - name of wavelet used for 1D transform
% order - order of the 2D wavelet transform (if biorthogonal -> order*10)
% orderz - order of 1D wavelet transform
% levels - number of levels in the wavelet transform
% tauMode - option for setting different wavelet thresholds (default 1)
% filtType - 'ht' or 'median'. 'ht' uses the standard algorithm with hard
%    thresh. on the first step, 'median' alternatively uses a median filter
% blockSize - size of the reference blocks in first step
% blockSizeWie - size of the reference blocks in the second step
% matchSize - size of search window for finding matching blocks
% cycleSpin - number of translations in the wavelet thresholding
% matchSpins - translation values for generating new sets of ref. blocks
% Wiener - true or false to perform the second step