function [U,Vavg,out] = denoiseAndAggregateGPU(I,S,sigma,opts)

% collaborative image denoising based on matched blocks
% inputs:
%   I - noisy grayscale image, values in [0,1]
%   S - cell matrix containing indices of matched blocks
%   d - block size
%   sigma - noise standard deviation
%   opts - options

% output: U - denoised image

% tau = 2.9*sigma; % hard thresholding constant
% for i = 1:opts.levels
%    % both of these seem to work fine
%    tau(i) = sigma*(3.3-.3*(i-1));
%    % tau(i) = 2.9*sigma;
% end



% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
[m,n] = size(I);
d = opts.blockSize;
Volume = zeros(m,n,opts.numMax); % initialize volume
Ipad = zeros(2*m,2*n); % make padded image
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard wavelet thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% wavelet decomposition
% opts.spin: number of spins in the wavelet cycle spinning
normSigma = sum(col(sigma));
if strcmp(opts.filtType,'ht')
    Vavg = zeros(m,n,opts.numMax); % volume that averages cycle spinning results
    for spin = 0:opts.cycleSpin-1
        % rotate volume and pass to GPU
        Vrec = gpuArray(single(circshift(Volume,[spin,spin,spin])));

        % get 3D wavelet filters and decompose
        Psi = getWaveFilters3D(opts.wname,opts.wnamez,opts.order,opts.orderz); 
        C = myWavDec3DFFT(Psi,opts.levels,Vrec);   

        % hard thresholding
        PhiP = gpuArray(single(ones(m,n,opts.numMax)));
        FPsi1 = fftn(Psi(:,:,:,1),[m,n,opts.numMax]);
        for i = 1:opts.levels 
            for j = 1:7
                waveFilter = sum(abs(fftn(Psi(:,:,:,j+1),[m,n,opts.numMax]).*PhiP),3);
                % [~,Zmax] = max(col(sum(waveFilter,3)));
                % tau = max(sigma(Zmax)*4,.1)
                % tau = 1;
                lambda = 1e2/normSigma;
                waveSigma = sum(col(waveFilter.*sigma).^2)/m/n;
                tau = lambda*sqrt(waveSigma);
                if j==3
                    tau = tau;% tau = .5;
                end
                C{i,j} = C{i,j}.*(abs(C{i,j})>tau);
            end
            PhiP = PhiP.*FPsi1;
        end

        % wavelet recontruction
        Vrec = myWavRec3DGPUfast(Psi,opts.levels,C);
        Vavg = Vavg + gather(circshift(Vrec/opts.cycleSpin,[-spin,-spin,-spin]));
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
    


% aggregate results, looping over each block
U = zeros(m+d,n+d); % pad for block matches that wrap around
out.K = zeros(m+d,n+d);
out.Z = zeros(m,n);
for y = 1:M
    for x = 1:N
        Nm = numel(S{y,x});
        % aggregate the blocks back into image and keep track of weights
        blocks = Vavg((y-1)*d+1:y*d,(x-1)*d+1:x*d,1:Nm);
        weight = 1/diff_norm_kZ(blocks,1,3,1);
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