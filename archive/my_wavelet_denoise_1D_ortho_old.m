function U = my_wavelet_denoise_1D_ortho(wname,levels,U,tau)

% Written by Toby Sanders
% orthonormal wavelet decomposition and reconstruction

% these are designed for 1-D problems

% inputs:
% wname - name of the wavelet used
% levels - number of levels in the decomposition
% U - signal to be denoised
% tau - thresholding value


[p,q,r] = size(U);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet decomposition step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = cell(levels+1,1); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    C{i} = circshift(imfilter(U,Hi_d','conv'),N/2);
    C{i} = C{i}(1:2:end);
    
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,Lo_d','conv'),N/2);
    U = U(1:2:end);
end
C{levels+1} = U; % last step returns approximation coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% throw away small detail coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:levels
    C{i} = max(abs(C{i})-tau,0).*sign(C{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet reconstruction algorithm, which performs the adjoint operations
% to the above loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(p,1); % store restored signal into U
Lo_d = flipud(Lo_d'); % flip filters since we need flipped conv.
Hi_d = flipud(Hi_d');

for i = 1:levels
   tmp = zeros(numel(C{i})*2,1);
   tmp(1:2:end) = C{i};
   tmp = circshift(imfilter(tmp,Hi_d,'conv'),-N/2+1);
   for j = 1:i-1
       tmp2 = zeros(numel(tmp)*2,1);
       tmp2(1:2:end) = tmp;
       tmp =  circshift(imfilter(tmp2,Lo_d,'conv'),-N/2+1);
   end
   U = U + tmp(1:p);
end

tmp = C{levels+1};
for i = 1:levels
    tmp2 = zeros(numel(tmp)*2,1);
    tmp2(1:2:end) = tmp;
    tmp = circshift(imfilter(tmp2,Lo_d,'conv'),-N/2+1);
end
U = U + tmp(1:p);


