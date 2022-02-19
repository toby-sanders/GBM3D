function [Psi,Psi2] = getWaveFilters3Dbior

% get 3D wavelet filters for GBM3D 
% biorthogonal order 1.5 in spatial dimensions crossed with Haar wavelets
% in time dimension


% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

% size of wavelet filters
N = 10;
Nz = 2;

% constants needed for wavelet filters
Ca = 0.016572815184060;
Cb = 0.121533978016438;
Cc = sqrt(2)/2; 

% bi-orthogonal decomposition wavelet filters
Lo_dxy = [Ca,-Ca,-Cb,Cb,Cc,Cc,Cb,-Cb,-Ca,Ca];
Hi_dxy = zeros(1,10); 
Hi_dxy(5) = -Cc; Hi_dxy(6) = Cc;

% Haar wavelet filters
Lo_d = [Cc,Cc];
Hi_d = [-Cc,Cc];


% build all 8 orthogonal 3D filters from the 1D filters
Psi = zeros(N,N,Nz,8,'single');
Psi2 = zeros(N,N,Nz,8,'single');
if gpuDeviceCount>0
    Psi = gpuArray(Psi);
    Psi2 = gpuArray(Psi2);
end

% decomposition filters
Psi(:,:,:,1) =  Lo_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,2) =  Hi_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,3) =  Lo_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,4) =  Lo_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,5) =  Hi_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,6) =  Hi_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,7) =  Lo_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,8) =  Hi_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);

% recontruction filters
Lo_dxy = abs(Hi_dxy);
Hi_dxy = [-Ca,-Ca,Cb,Cb,-Cc,Cc,-Cb,-Cb,Ca,Ca];
Psi2(:,:,:,1) =  Lo_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi2(:,:,:,2) =  Hi_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi2(:,:,:,3) =  Lo_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi2(:,:,:,4) =  Lo_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi2(:,:,:,5) =  Hi_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi2(:,:,:,6) =  Hi_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi2(:,:,:,7) =  Lo_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);
Psi2(:,:,:,8) =  Hi_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);


