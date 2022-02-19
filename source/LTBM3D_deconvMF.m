function [U,out] = LTBM3D_deconv(I,hhat,sigma,opts)

sZ = 32;
LTopts.blockSize = sZ;
LTopts.numMax = 16; % max number of blocks to match for each ref. block
LTopts.numMin = 16; % min number of blocks to match for each ref. block
LTopts.wname = 'bior';
LTopts.wnamez = 'db';
LTopts.order =13;
LTopts.orderz = 1;
LTopts.levels = 3;
LTopts.cycleSpin = 2;
LTopts.matchSpins =  [0 4 sZ/2-1];
LTopts.wiener = false;
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

filt = sum(conj(hhat),3)./(sum(hhat2,3) + lambda*V);
out.recWei2 = real(ifft2(sum(FI.*conj(hhat),3)./(sum(hhat2,3) + lambda*V)));
out.sigmaPSD = sigma^2*abs(filt).^2;
U = LTBM3D(out.recWei2,out.sigmaPSD,LTopts);
