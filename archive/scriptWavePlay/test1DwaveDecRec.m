clear;
d = 128;
levels = 3;

load('bior15Filters');
Lo_dxy = LoD;
Hi_dxy = HiD;
N = numel(LoD);


Psi = [-sqrt(2)/2;sqrt(2)/2]; Psi2 = Psi;
Phi = abs(Psi); Phi2 = Phi;


load('bior15Filters');
Phi = LoD';
Psi = HiD';
Psi2 = flipud(HiR');
Phi2 = flipud(LoR');



U = zeros(d,1);
U(d/2) = 1;

C = myWavDec1DFFT(Psi,Phi,levels,U);

U2 = myWavRec1DGPUfast(Psi2,Phi2,levels,C);

s = 0;
for i = 1:numel(C)
    s = s + norm(C{i}(:))^2;
end
s
norm(U(:))