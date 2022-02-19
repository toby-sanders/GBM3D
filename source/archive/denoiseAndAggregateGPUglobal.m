function [U,Vavg,out] = denoiseAndAggregateGPU(I,S,sigma,opts)
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

% tau = 2.9*sigma; % hard thresholding constant
% for j = 1:7
% for i = 1:opts.levels
   % both of these seem to work fine
  %  tau(i,j) = sigma*(3.3-.3*(i-1) + j*.1);
  % tau(i,j) = 2.4*sigma;
% end
% end
tau = sigma*getTau(opts.levels,opts.tauMode);

% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
[m,n] = size(I);
d = opts.blockSize;
Volume = zeros(m,n,opts.numMax,'single'); % initialize volume
Ipad = zeros(2*m,2*n,'single'); % make padded image
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
for i = 1:M % loop over each block
    for j = 1:N
        xR = Ipad((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        for k = 1:numel(S{i,j})
            [a,b] = ind2sub([m,n],S{i,j}(k));
            Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,k) = Ipad(a:a+d-1,b:b+d-1);
        end
        % pad the rest of the block in the volume with the reference block
        Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,k+1:opts.numMax) = repmat(xR,1,1,opts.numMax-k);
    end
end

% output this just for prototyping
out.dd = squeeze(sum(sum(abs(diff(Volume,1,3)))))./...
    squeeze(sum(sum(abs(Volume(:,:,1:end-1)))));



% initialize GPU variabes and wavelet filters
Vrec = zeros(m,n,opts.numMax,opts.cycleSpin,'gpuArray');
Volume = gpuArray(Volume);
if strcmp(opts.wname,'bior')
    [Psi,Psi2] = getWaveFilters3Dbior(opts.wnamez,opts.order,opts.orderz);
else
    Psi = getWaveFilters3D(opts.wname,opts.wnamez,opts.order,opts.orderz);
    Psi2 = Psi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard wavelet thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% opts.cycleSpin: number of spins in the wavelet cycle spinning
if strcmp(opts.filtType,'ht')
   % Vavg = zeros(m,n,opts.numMax); % volume that averages cycle spinning results
    for spin = 0:opts.cycleSpin-1
        % rotate volume for cyclespinning
        Vrec(:,:,:,spin+1) = circshift(Volume,[spin,spin,spin]);

        % wavelet decomposition
        C = myWavDec3DFFT(Psi,opts.levels,Vrec(:,:,:,spin+1));   

        % hard thresholding
        for i = 1:opts.levels 
            for j = 1:7
                C{i,j} = C{i,j}.*(abs(C{i,j})>tau(i,j));
            end
        end

        % wavelet recontruction
        Vrec(:,:,:,spin+1) = circshift(myWavRec3DGPUfast(Psi2,opts.levels,C),[-spin,-spin,-spin]);
        % Vavg = Vavg + gather(circshift(Vrec/opts.cycleSpin,[-spin,-spin,-spin]));
    end
else % Wiener filter on DCT of volume
    Volume = dct(Volume,[],1);
    Volume = dct(Volume,[],2);
    Volume = dct(Volume,[],3);
    filt = abs(Volume).^2./(abs(Volume).^2 + sigma^2);
    Volume = Volume.*filt;
    Volume = idct(Volume,[],3);
    Volume = idct(Volume,[],2);
    Vavg = idct(Volume,[],1);
end
Vavg = gather(mean(Vrec,4));


% aggregate results, looping over each block
U = zeros(m+d,n+d); % pad for block matches that wrap around
out.K = zeros(m+d,n+d);
out.Z = zeros(m,n);
D = getBM3DAggWeights(3,1);
for y = 1:M
    for x = 1:N
        Nm = numel(S{y,x});
        % aggregate the blocks back into image and keep track of weights
        blocks = Vavg((y-1)*d+1:y*d,(x-1)*d+1:x*d,1:Nm);
        weight = 1/D(blocks);
        % weight = 1/diff_norm_kZ(blocks,1,3,1);
        % weight = 1;
        for i = 1:Nm
            [a,b] = ind2sub([m,n],S{y,x}(i));
            U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight*blocks(:,:,i);
            out.K(a:a+d-1,b:b+d-1) = out.K(a:a+d-1,b:b+d-1) + weight;
        end
        out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = Nm;
    end
end

% reweight image by number of times each pixel is visited
U = U./out.K;
U = U(1:m,1:n); % delete the padding
out.K = out.K(1:m,1:n);