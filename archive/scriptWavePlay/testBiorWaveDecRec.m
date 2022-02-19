clear;
levels = 2;
d0 = 256;
dim = 1;

% load(['bior15Filters']);
% Phi = LoD';
% Psi = HiD';
% Psi2 = (Psi);
% Phi2 = (Phi);
load('bior15Filters');
Phi = LoD';
Psi = HiD';
Psi2 = flipud(HiR');
Phi2 = flipud(LoR');

% Psi = [-sqrt(2)/2;sqrt(2)/2]; Psi2 = Psi;
% Phi = abs(Psi); Phi2 = Phi;
x = zeros(d0,1);
x(d0/2) = 1;
x0 = x;
U = x;

C = zeros(d0,1);
cnt = 0;
for i = 1:levels
        % high pass filter and downsample
    FU = fft(U,size(U,dim),dim);
    d = ifft(FU.*fft(Psi,size(U,dim),dim),size(U,dim),dim);
    C(cnt+1:cnt+size(d,dim)/2) = real(d(1:2:end));
    cnt = cnt + size(d,dim)/2;
    
    % low pass filter, downsample, and move to next level
    % U = circshift(imfilter(U,Lo_d','conv','circular'),N/2);
    U = real(ifft(FU.*fft(Phi,size(U,dim),dim),size(U,dim),dim));
    U = U(1:2:end);
end
C(cnt+1:end) = U; % last step returns approximation coefficients




Cup = cell(levels+1,1);
cnt = 0;
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = zeros(d0/2^(i-1),1);
    Cup{i}(1:2:end) = C(cnt+1:cnt+d0/2^i);
    cnt = cnt + d0/2^i;
    FPsi = conj(fft(Psi2,size(Cup{i},dim),dim));
    Cup{i} = real(ifft(fft(Cup{i},size(Cup{i},dim),dim).*FPsi,size(Cup{i},dim),dim));
end
Cup{levels+1} = C(cnt+1:end);

% iteratively LPF and upsample
for i = levels+1:-1:2
    tmp = zeros(size(Cup{i-1},dim),1,'single');
    tmp(1:2:end) = Cup{i};
    FPsi = conj(fft(Phi2,size(tmp,dim),dim));
    Cup{i-1} = Cup{i-1} + real(ifft(fft(tmp,size(tmp,dim),dim).*FPsi,size(tmp,dim),dim));
end
U = Cup{1}; % next level coefficients added


figure(111);plot(U);