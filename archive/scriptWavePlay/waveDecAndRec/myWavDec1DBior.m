function C = myWavDec1DBior(LoD,HiD,levels,U)

% orthogonal wavelet transform
% i.e. wavelet decomposition algorithm
N = numel(LoD);

C = cell(levels+1,1); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    d = circshift(imfilter(U,HiD','conv','circular'),N/2); % detail coef.
    C{i} = d(1:2:end);
    
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,LoD','conv','circular'),N/2);
    U = U(1:2:end);
end
C{levels+1} = U; % last step returns approximation coefficients
