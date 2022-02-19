function opts = setBM3Dopts(profile)

% set GBM3D parameters based on the selected profile
% options for the profile are:
%   - 'default' (accurate and somewhat fast)
%   - 'fast'    (fast)
%   - 'superFast'   (fastest)
%   - 'accuracy'  (most accurate)

% Written by Toby Sanders @Lickenbrock Tech
% 3-17-2021
if ~sum(strcmp(profile,{'default','fast','superFast','accuracy'}))
    profile = 'default';
end

% set parameters
switch profile
    case 'default'
        opts.numMax = 16;
        opts.numMin = 16;
        opts.wname = 'bior';
        opts.wnamez = 'db';
        opts.order = 15;
        opts.orderz = 1;
        opts.levels = 3;
        opts.cycleSpin = 2;
        opts.matchSpins =  [0 4]; % two spins on both wavelet and Wiener
        opts.matchSpinsWie = [0 4];
        opts.tauMode = 1;
        opts.blockSize = 16;
        opts.blockSizeWie = 8;
        opts.matchSize = 32;
        opts.Wiener = true;
        opts.filtType = 'ht';
    case 'fast'
        opts.numMax = 16;
        opts.numMin = 16;
        opts.wname = 'bior';
        opts.wnamez = 'db';
        opts.order = 15;
        opts.orderz = 1;
        opts.levels = 3;
        opts.cycleSpin = 2;
        opts.matchSpins =  [0 4 6]; % add one spin
        opts.matchSpinsWie = [0]; % no spins on Wiener step
        opts.tauMode = 1;
        opts.blockSize = 16;
        opts.blockSizeWie = 8;
        opts.matchSize = 32;
        opts.Wiener = true; % Perform Wiener step
        opts.filtType = 'ht';
    case 'superFast'
        opts.numMax = 16;
        opts.numMin = 16;
        opts.wname = 'bior';
        opts.wnamez = 'db';
        opts.order = 15;
        opts.orderz = 1;
        opts.levels = 3;
        opts.cycleSpin = 2;
        opts.matchSpins =  [0 4 6]; % one less spin
        opts.matchSpinsWie = [0];
        opts.tauMode = 1;
        opts.blockSize = 32; % larger block Size
        opts.blockSizeWie = 8;
        opts.matchSize = 32;
        opts.Wiener = false; % skip Wiener step
        opts.filtType = 'ht';
    case 'accuracy'
        opts.numMax = 16;
        opts.numMin = 16;
        opts.wname = 'bior';
        opts.wnamez = 'db';
        opts.order = 15;
        opts.orderz = 1;
        opts.levels = 3;
        opts.cycleSpin = 2;
        opts.matchSpins =  [0 4 6]; % three spins on both wavelet and Wiener
        opts.matchSpinsWie = [0 4 6];
        opts.tauMode = 1;
        opts.blockSize = 16;
        opts.blockSizeWie = 8;
        opts.matchSize = 32;
        opts.Wiener = true;
        opts.filtType = 'ht';
end
        
        
        