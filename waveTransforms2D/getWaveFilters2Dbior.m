function [Psi,Psi2] = getWaveFilters2Dbior(order)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

fname = ['bior',num2str(order),'Filters'];
load(['waveletFilters',filesep,fname]);
% load(fname);
N = numel(LoD);
Lo_dxy = LoD;
Hi_dxy = HiD;

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


