function C = myWavDec(wname,levels,U)


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
