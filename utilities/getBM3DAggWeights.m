function D = getBM3DAggWeights(kx,kz)



% finite difference operators for higher order TV
% k is the order of the transform
%
% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

D = @(U)D_Forward(U,kx,kz);


% high order finite differences
function W = D_Forward(U,kx,kz)

W = sum(col(abs(diff(U,kx,2))))*2^(1-kx);
W = W + sum(col(abs(diff(U,kx,1))))*2^(1-kx);
W = W + sum(col(abs(diff(U,kx,3))))*2^(1-kz);

    
    