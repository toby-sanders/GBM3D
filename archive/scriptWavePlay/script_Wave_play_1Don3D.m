clear;
d = 128;
SNR = 3;
levels = 3;
order = 1;
tau = 0.0;
rng(2021);

addpath('waveDecAndRec');
wname = ['db',num2str(order)];

P = phantom(d);
P = repmat(P,1,1,d);
P = add_Wnoise(P,SNR);
% x = sin(2*pi*linspace(0,1,d));
% P = reshape(x,1,1,d);

% C1 = myWavDec1D(wname,levels,P(:));
% P1 = myWavRec1DFast(wname,levels,C1);

C = myWavDec1Don3D(wname,levels,P);
P2 = myWavRec1Don3D(wname,levels,C);


myrel(P2,P)