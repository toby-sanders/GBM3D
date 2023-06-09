function [U,out] = Tikhonov_SUPER(I,hhat,V)

% This function solves
% a tikhonov regularized super resolution deconvolution problem

% INPUTS:
% "I" is the low res. image
% hhat is the Fourier transform of the PSF (at high res.), "coherent
% transfer function"
% K is the upsampling factor
% V is the regularization operator (in Fourier domain)

% The solver uses conjugate gradient, and it is solved in Fourier domain,
% which reduces the number of FFTs needed in each iteration from 4 to 2


% Written by Toby Sanders @Lickenbrock Tech.
% 01/19/2022

% set image dimensions
K = 2;
[m,n] = size(I);
M = m*K;
N = n*K;


% opts = check_tik_opts(opts);
tol = 1e-4;
iter = 50;

g = zeros(M,N);
g([1:K],[1:K]) = 1/K^2;
g = fraccircshift(g,[-K/2 + 1/2, -K/2 + 1/2]);
ghat = fft2(g);
ghat = ghat.*hhat;


if ~exist('V','var')
    V = my_Fourier_filters(opts.order,1,M,N,1);
end

blockA = abs(ghat(:,1:N/2)).^2;
blockD = abs(ghat(:,N/2+1:end)).^2;
blockB = conj(ghat(:,1:N/2)).*ghat(:,N/2+1:end);



