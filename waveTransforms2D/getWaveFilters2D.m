function Psi = getWaveFilters2D(wname,order)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


fname = [wname,num2str(order),'Filters'];
load([fname]);
N = numel(Lo_d);

% build all 8 orthogonal 3D filters from the 1D filters
Psi = gpuArray(zeros(N,N,4,'single'));
Psi(:,:,1) =  Lo_d'*Lo_d;
Psi(:,:,2) =  Hi_d'*Lo_d;
Psi(:,:,3) =  Lo_d'*Hi_d;
Psi(:,:,4) =  Hi_d'*Hi_d;




