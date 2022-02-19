function U = myWavRec(wname,levels,C)


[p,q,r] = size(C);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet reconstruction algorithm, which performs the adjoint operations
% to the above loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(p,1); % store restored signal into U
Lo_d = Lo_d'; % flip filters since we need transpose operation
Hi_d = Hi_d';

cnt = 0;
for i = 1:levels
    cnt = cnt+p/(2^i);
end

% upsample the detail coefficients
Cup = zeros(cnt + p ,1);
Cup(1:2:2*cnt) = C(1:cnt);
Cup(2*cnt+1:end) = C(cnt+1:end);

% HPF the upsampled detail coefficients
Cup(1:2*cnt) = circshift(imfilter(Cup(1:2*cnt),Hi_d,'corr','replicate'),-N/2+1);
U = U + Cup(1:p); % first level coefficients added
Cup = Cup(p+1:end); % throw out these coefficients and proceed
for i = 1:levels
   % upsample coefficients and LPF
   tmp = zeros(numel(Cup)*2,1);
   tmp(1:2:end) = Cup;
   Cup =  circshift(imfilter(tmp,Lo_d,'corr','replicate'),-N/2+1);
   U = U + Cup(1:p); % next level coefficients added
   Cup = Cup(p+1:end); % throw out these coefficients and proceed
end
