function g = getSuperResGrad(hhat,mask,uhat,V)

% gradient operation need for CG-based Tikhonov solver in super resolution
% imaging

% written by Toby Sanders
% @ Lickenbrock Tech., Inc.
% 1/18/22

uhat = reshape(uhat,size(mask,1),size(mask,2));
g = V.*uhat; % no lambda needed, already embedded into V
g2 = ifft2(hhat.*uhat);
g2 = fft2(g2.*mask).*conj(hhat);
g = col(g + g2);