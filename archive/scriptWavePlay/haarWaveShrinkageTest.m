% For the second step of BM3D, I wanted to know if some of the coefficients
% in the 3D transform should NOT be filtered, i.e. the approximation
% coefficients

d = 16;
d0 = 512;
sigma = .1;

Psi = single(reshape([-sqrt(2)/2,sqrt(2)/2],1,1,2));
Phi = abs(Psi);

P = phantom(d0);
x = P(191:190+d,206:205+d);

y = repmat(x,1,1,d);
y0 = y;
y = y + randn(size(y))*sigma;


C = myWavDec1D(Psi,Phi,2,y,3);
C2 = dct(C,[],1);
C3 = dct(C,[],2);