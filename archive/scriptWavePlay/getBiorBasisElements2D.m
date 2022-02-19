clear;
levels = 2;
d0 = 256;


load('bior15Filters');
Lo_dxy = LoD;
Hi_dxy = HiD;
N = numel(LoD);

% build all 8 orthogonal 3D filters from the 1D filters
Psi = zeros(N,N,4,'single');
Psi2 = zeros(N,N,4,'single');

Psi(:,:,1) =  Lo_dxy'*Lo_dxy;
Psi(:,:,2) =  Hi_dxy'*Lo_dxy;
Psi(:,:,3) =  Lo_dxy'*Hi_dxy;
Psi(:,:,4) =  Hi_dxy'*Hi_dxy;


Lo_dxy = fliplr(LoR);
Hi_dxy = fliplr(HiR);

Psi2(:,:,1) =  Lo_dxy'*Lo_dxy;
Psi2(:,:,2) =  Hi_dxy'*Lo_dxy;
Psi2(:,:,3) =  Lo_dxy'*Hi_dxy;
Psi2(:,:,4) =  Hi_dxy'*Hi_dxy;


P = phantom(d0);

C = myWavDec2DFFT(Psi,levels,P);

P2 = myWavRec2DGPUfast(Psi2,levels,C);

myrel(P,P2)

s = 0;
for i = 1:numel(C)
    s = s + norm(C{i}(:))^2;
end
s
norm(P(:))^2