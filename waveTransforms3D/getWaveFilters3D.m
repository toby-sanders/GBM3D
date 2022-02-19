function Psi = getWaveFilters3D(wname,wnamez,order,orderz)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

fname = [wname,num2str(order),'Filters'];
load([fname]);
% load(fname);
N = numel(Lo_d);
Lo_dxy = Lo_d;
Hi_dxy = Hi_d;

fname = [wnamez,num2str(orderz),'Filters'];
load(['waveletFilters',filesep,fname]);
% load(fname);
Nz = numel(Lo_d);

% build all 8 orthogonal 3D filters from the 1D filters
Psi = zeros(N,N,Nz,8,'single');
if gpuDeviceCount>0
    Psi = gpuArray(Psi);
end

Psi(:,:,:,1) =  Lo_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,2) =  Hi_dxy'*Lo_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,3) =  Lo_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,4) =  Lo_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,5) =  Hi_dxy'*Hi_dxy.*reshape(Lo_d,1,1,Nz);
Psi(:,:,:,6) =  Hi_dxy'*Lo_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,7) =  Lo_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);
Psi(:,:,:,8) =  Hi_dxy'*Hi_dxy.*reshape(Hi_d,1,1,Nz);


