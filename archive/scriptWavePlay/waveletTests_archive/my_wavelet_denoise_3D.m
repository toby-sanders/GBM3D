function U = my_wavelet_denoise_3D(wname,levels,U)

% Written by Toby Sanders

% these are designed for 3-D problems
[p,q,r] = size(U);
% [Lo_d,Hi_d] = wfilters(wname);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
lo_x = fft(Lo_d,q);
lo_y = fft(Lo_d,p);
lo_z = fft(Lo_d,r);
hi_x = fft(Hi_d,q);
hi_y = fft(Hi_d,p);
hi_z = fft(Hi_d,r);

[LOX,LOY,LOZ] = meshgrid(lo_x,lo_y,lo_z);
[HIX,HIY,HIZ] = meshgrid(hi_x,hi_y,hi_z);


U = reshape(U,p,q,r);
Ux = fft(U,q,2);
Uy = fft(U,p,1);
Uz = fft(U,r,3);

cx = zeros(p,q,r,levels+1);
cy = zeros(p,q,r,levels+1);
cz = zeros(p,q,r,levels+1);
for i = 1:levels
    cx(:,:,:,i) = ifft(Ux.*HIX,q,2);
    Ux = Ux.*LOX;
    
    cy(:,:,:,i) = ifft(Uy.*HIY,p,1);
    Uy = Uy.*LOY;
    
    cz(:,:,:,i) = ifft(Uz.*HIZ,r,3);
    Uz = Uz.*LOZ;
end
cx(:,:,:,levels+1) = ifft(Ux,q,2);
cy(:,:,:,levels+1) = ifft(Uy,p,1);
cz(:,:,:,levels+1) = ifft(Uz,r,3);

tau = 000;
cx(:,:,:,1:levels) = max(abs(cx(:,:,:,1:levels))-tau,0).*sign(cx(:,:,:,1:levels));
cy(:,:,:,1:levels) = max(abs(cy(:,:,:,1:levels))-tau,0).*sign(cy(:,:,:,1:levels));
cz(:,:,:,1:levels) = max(abs(cz(:,:,:,1:levels))-tau,0).*sign(cz(:,:,:,1:levels));

cx = fft(cx,q,2);
cy = fft(cy,p,1);
cz = fft(cz,r,3);
Ux = zeros(p,q,r);
Uy = zeros(p,q,r);
Uz = zeros(p,q,r);

LHFx = HIX;
LHFy = HIY;
LHFz = HIZ;
for i = 1:levels
    Ux = Ux + ifft(cx(:,:,:,i).*conj(LHFx),q,2)/(2^i);
    Uy = Uy + ifft(cy(:,:,:,i).*conj(LHFy),p,1)/(2^i);
    Uz = Uz + ifft(cz(:,:,:,i).*conj(LHFz),r,3)/(2^i);
    LHFx = LHFx.*LOX;
    LHFy = LHFy.*LOY;
    LHFz = LHFz.*LOZ;
    
end
Ux = Ux + ifft(cx(:,:,:,levels+1).*conj(LOX.^levels),q,2)/(2^levels);
Uy = Uy + ifft(cy(:,:,:,levels+1).*conj(LOY.^levels),p,1)/(2^levels);
% Uz = Uz + ifft(cz(:,:,:,levels+1).*conj(LOZ.^levels),r,3)/(2^levels);
U = (Ux + Uy)/2;% + Uz)/3;



