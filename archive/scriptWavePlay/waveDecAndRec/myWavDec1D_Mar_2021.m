function C = myWavDec(wname,levels,U)

% orthogonal wavelet transform
% i.e. wavelet decomposition algorithm
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

C = cell(levels+1,1); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    FU = fft(U);
    d = ifft(FU.*fft(Hi_d',size(U,1)));
    % d = circshift(imfilter(U,Hi_d','conv','circular'),N/2); % detail coef.
    C{i} = real(d(1:2:end));
    
    % low pass filter, downsample, and move to next level
    % U = circshift(imfilter(U,Lo_d','conv','circular'),N/2);
    U = real(ifft(FU.*fft(Lo_d',size(U,1))));
    U = U(1:2:end);
end
C{levels+1} = U; % last step returns approximation coefficients
