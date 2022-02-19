function [U,Vavg,out] = denoiseAndAggregateGPUColored(Volume,I,S,sigmaPSD,opts)

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

tau = getColoredTau(opts.levels,opts.order,sigmaPSD);

% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
[m,n] = size(I);
d = opts.blockSize;

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
Vavg = min(max(gather(mean(Vrec,4)),0),1);


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
            [a,b] = ind2sub([d,d],S{y,x}(i));
            a = a + (y-3/2)*d;
            b = b + (x-3/2)*d;
            if a>0 & b>0
            U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight*blocks(:,:,i);
            out.K(a:a+d-1,b:b+d-1) = out.K(a:a+d-1,b:b+d-1) + weight;
            end
        end
        out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = Nm;
    end
end

% reweight image by number of times each pixel is visited
U = U./out.K;
U = U(1:m,1:n); % delete the padding
out.K = out.K(1:m,1:n);