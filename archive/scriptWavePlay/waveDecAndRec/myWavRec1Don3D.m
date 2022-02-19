function U = myWavRec1D(wname,levels,C)

% orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
[p,q,r] = size(C{1});
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

U = zeros(p,q,2*r); % store restored signal into U
Psi = reshape(Hi_d,1,1,numel(Hi_d));
Phi = reshape(Lo_d,1,1,numel(Lo_d));


Cup = cell(levels+1,1);
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = zeros(p,q,2*size(C{i},3));
    Cup{i}(:,:,1:2:end) = C{i};
    FPsi = conj(fft(Psi,size(Cup{i},3),3));
    Cup{i} = real(ifft(fft(Cup{i},size(Cup{i},3),3).*FPsi,size(Cup{i},3),3));
end
Cup{levels+1} = C{levels+1};

% iteratively LPF and upsample
for i = levels+1:-1:2
    tmp = zeros(p,q,size(Cup{i-1},3),'single');
    tmp(:,:,1:2:end) = Cup{i};
    FPsi = conj(fft(Phi,size(tmp,3),3));
    Cup{i-1} = Cup{i-1} + real(ifft(fft(tmp,size(tmp,3),3).*FPsi,size(tmp,3),3));
end
   U = Cup{1}; % next level coefficients added