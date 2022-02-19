function [U,out] = Tikhonov_SUPER(I,K,hhat,V)

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

c = zeros(M,N);
c(1:K:end,1:K:end) = I;
c = conj(ghat).*fft2(c);
mask = zeros(M,N);
mask(1:K:end,1:K:end) = 1;
B = @(x)getSuperResGrad(ghat,mask,x,V);

% A(b,2);
tic;
[U,out] = my_local_cgs(B,c(:),iter,tol);
out.total_time = toc;
U = real(ifft2(reshape(U,M,N)));

function [x,out] = my_local_cgs(A,b,iter,tol)
% CG algorithm from wiki

x = zeros(size(b));
out = [];
r = b-A(x);
p = r;
rel_chg = zeros(iter,1);
for i = 1:iter
    Ap = A(p);
    alpha = r'*r./(p'*Ap);
    xp = x;
    x = x + alpha*p;
    rp = r;
    r = r-alpha.*Ap;
    rel_chg(i) = real(alpha^2*(p'*p)/(x'*x));
    if rel_chg(i) < tol^2; break;end
    beta = r'*r./(rp'*rp);
    p = r+beta.*p;
end
out.iters = i;
out.rel_chg = sqrt(rel_chg(1:i));
