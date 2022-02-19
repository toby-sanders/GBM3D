function V = my3DWiener(V,p)

if nargin<2, p = 1; end
[x,y,z] = size(V);

f1 = (cos(2*pi*linspace(0,1,x))+1)/2;f1 = f1';
f2 = (cos(2*pi*linspace(0,1,y))+1)/2;
f3 = reshape((cos(2*pi*linspace(0,1,z))+1)/2,1,1,z);

F = ones(x,y,z);
F = F.*(f1*f2);
F = F.*f3;
V = real(ifft2(F.^p.*fft2(V)));