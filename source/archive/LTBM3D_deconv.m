function [U,out] = LTBM3D_deconv(I,hhat,sigma,opts)

sZ = 32;
LTopts.numMax = 16;
LTopts.numMin = 16;
LTopts.wname = 'bior';
LTopts.wnamez = 'db';
LTopts.order = 15;
LTopts.orderz = 1;
LTopts.levels = 3;
LTopts.cycleSpin = 2;
LTopts.matchSpins =  [0 4];
LTopts.tauMode = 1;
LTopts.blockSize = 16;
LTopts.blockSizeWie = 8;
LTopts.matchSize = 32;
LTopts.Wiener = true;
LTopts.filtType = 'ht';


if ~isfield(opts,'order'),opts.order = 1; end
if ~isfield(opts,'levels'), opts.levels = 1; end
if ~isfield(opts,'lambda'), opts.lambda = 200; end


[d1,d2,nF] = size(I);
if size(hhat,1)~=d1 || size(hhat,2)~=d2
    error('hhat should have same dimensions at I');
end

FI = fft2(I);
hhat2 = abs(hhat).^2;

V = my_Fourier_filters(opts.order,opts.levels,d1,d2,1);
parm = opts;
parm.theta = sigma^2*opts.lambda;
[out.recWei,out.ME] = HOTVL2_deblurMF(ifft2(hhat),I,parm);
lambda = out.ME.thetas(end)/2

filt = conj(hhat)./(hhat2 + lambda*V);
out.recWei2 = real(ifft2(filt.*FI));
out.sigmaPSD = sigma^2*abs(filt).^2;
U = GBM3D(out.recWei2,out.sigmaPSD,LTopts);
