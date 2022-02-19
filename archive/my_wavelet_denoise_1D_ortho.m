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
C = zeros(p,1); % coef. at each level
cnt = 0;
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    d = circshift(imfilter(U,Hi_d','conv','circular'),N/2); % detail coef.
    C(1+cnt:cnt+p/(2^i)) = d(1:2:end);
    cnt = cnt+p/(2^i);
    
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,Lo_d','conv','circular'),N/2);
    U = U(1:2:end);
end
C(cnt+1:end) = U; % last step returns approximation coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% throw away small detail coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C(1:cnt) = max(abs(C(1:cnt))-tau,0).*sign(C(1:cnt));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet reconstruction algorithm, which performs the adjoint operations
% to the above loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(p,1); % store restored signal into U
Lo_d = Lo_d'; % flip filters since we need transpose operation
Hi_d = Hi_d';


% upsample the detail coefficients
Cup = zeros(cnt + p ,1);
Cup(1:2:2*cnt) = C(1:cnt);
Cup(2*cnt+1:end) = C(cnt+1:end);

% HPF the upsampled detail coefficients
Cup(1:2*cnt) = circshift(imfilter(Cup(1:2*cnt),Hi_d,'corr','circular'),-N/2+1);
U = U + Cup(1:p); % first level coefficients added
Cup = Cup(p+1:end); % throw out these coefficients and proceed
for i = 1:levels
   % upsample coefficients and LPF
   tmp = zeros(numel(Cup)*2,1);
   tmp(1:2:end) = Cup;
   Cup =  circshift(imfilter(tmp,Lo_d,'corr','circular'),-N/2+1);
   U = U + Cup(1:p); % next level coefficients added
   Cup = Cup(p+1:end); % throw out these coefficients and proceed
end


