clear;
levels = 3;
d0 = 64;
dim = 1;


load('bior15Filters');
Phi = LoD';
Psi = HiD';
Psi2 = flipud(HiR');
Phi2 = flipud(LoR');

% Psi = [-sqrt(2)/2;sqrt(2)/2]; Psi2 = Psi;
% Phi = abs(Psi); Phi2 = Phi;
x = zeros(d0,1);
x(1) = 1;
x0 = x;
M = zeros(d0);
for j = 1:d0

U = zeros(d0,1);
U(j) = 1;



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

M(:,j) = C;
end


M2 = M;
for j = 1:d0
    
    C = zeros(d0,1);
    C(j) = 1;
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
M2(:,j) = U;
end

Tforward =[ 0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110
               -0.225454819240296  -0.461645582253923  -0.461645582253923  -0.225454819240296   0.225454819240296   0.461645582253923   0.461645582253923   0.225454819240296
                0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840  -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280
               -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280   0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840
                0.707106781186550  -0.707106781186550                   0                   0                   0                   0                   0                   0
                                0                   0   0.707106781186550  -0.707106781186550                   0                   0                   0                   0
                                0                   0                   0                   0   0.707106781186550  -0.707106781186550                   0                   0
                                0                   0                   0                   0                   0                   0   0.707106781186550  -0.707106781186550];
                            
                            
figure(177);
subplot(2,2,1);imagesc(M);colorbar;
subplot(2,2,2);imagesc(Tforward);colorbar;
subplot(2,2,3);imagesc(M2);colorbar;
subplot(2,2,4);imagesc(inv(M));colorbar;


D = zeros(size(M));
for i = 1:size(D,1)
    D(i,i) = randn(1);
    if D(i,i)<0
        D(i,i) = 0;
    else
        D(i,i) = 1;
    end
end

T = M2*D*M;
figure(178);
subplot(1,2,1);imagesc(abs(T));
subplot(1,2,2);imagesc(abs(T'));



myrel(T,T')
myrel(M2,inv(M))