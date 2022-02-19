function w = makeInterpWindow(a,b,k)

% make interpolation windows for distributed version of GBM3D

% Written by Toby Sanders @Lickenbrock Tech
% 3-17-2021

w = ones(a,b);
ramp = linspace(.01,1,k);
w(1:k,:) = w(1:k,:).*repmat(ramp',1,b);
w(:,1:k) = w(:,1:k).*repmat(ramp,a,1);
w(end-k+1:end,:) = w(end-k+1:end,:).*flipud(repmat(ramp',1,b));
w(:,end-k+1:end) = w(:,end-k+1:end).*fliplr(repmat(ramp,a,1));