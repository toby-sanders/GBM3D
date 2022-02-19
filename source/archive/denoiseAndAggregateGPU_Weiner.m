function [U,Vavg,out] = denoiseAndAggregateGPU(I,U,S,sigma,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% collaborative image denoising based on matched blocks
% inputs:
%   I - noisy grayscale image, values in [0,1]
%   S - cell matrix containing indices of matched blocks
%   d - block size
%   sigma - noise standard deviation
%   opts - options

% output: U - denoised image

tau = sigma*getTau(opts.levels,opts.tauMode);

% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
[m,n] = size(I);
d = opts.blockSize;
Ipad = zeros(2*m,2*n,'single'); % make padded image
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
Ipad = Ipad(1:m+d,1:n+d);
Ipad = circshift(Ipad,[d/2,d/2]);
Upad = zeros(2*m,2*n,'single'); % make padded image
Upad(1:m,1:n) = U; Upad(1:m,n+1:end) = U; 
Upad(m+1:end,1:n) = I;Upad(m+1:end,n+1:end) = U;
Upad = Upad(1:m+d,1:n+d);
Upad = circshift(Upad,[d/2,d/2]);


% aggregate results, looping over each block
U = zeros(m+d,n+d); % pad for block matches that wrap around
out.K = zeros(m+d,n+d);
out.Z = zeros(m,n);
D = getBM3DAggWeights(3,1);
for i = 1:M % loop over each block
    for j = 1:N
        Nm = numel(S{i,j});
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        Itmp = Ipad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
        Utmp = Upad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
        V1 = zeros(d,d,Nm);V2 = V1;
        for k = 1:Nm
            [a,b] = ind2sub([d,d],S{i,j}(k));
            V1(:,:,k) = Itmp(a:a+d-1,b:b+d-1);
            V2(:,:,k) = Utmp(a:a+d-1,b:b+d-1);
        end
        V1 = dct(V1,[],1);
        V1 = dct(V1,[],2);
        V1 = dct(V1,[],3);
        V2 = dct(V2,[],1);
        V2 = dct(V2,[],2);
        V2 = dct(V2,[],3);
        filt = abs(V2).^2./(abs(V2).^2 + sigma^2);
        V = V1.*filt;
        V = idct(V,[],3);
        V = idct(V,[],2);
        blocks = idct(V,[],1);
        
        weight = 1/D(blocks);
        % weight = 1/diff_norm_kZ(blocks,1,3,1);
        % weight = 1;
        for k = 1:Nm
            [a,b] = ind2sub([d,d],S{i,j}(k));
            a = a + (i-3/2)*d;
            b = b + (j-3/2)*d;
            if a>0 & b>0
            U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight*blocks(:,:,k);
            out.K(a:a+d-1,b:b+d-1) = out.K(a:a+d-1,b:b+d-1) + weight;
            end
        end
        out.Z((i-1)*d+1:i*d,(j-1)*d+1:j*d) = Nm;
    end
end


% reweight image by number of times each pixel is visited
U = U./out.K;
U = U(1:m,1:n); % delete the padding
out.K = out.K(1:m,1:n);