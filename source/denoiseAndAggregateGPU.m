function [U,W] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% collaborative image denoising based on matched blocks


% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
d = opts.blockSize;
m = M*d; n = N*d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard wavelet thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% opts.cycleSpin: number of spins in the wavelet cycle spinning
if strcmp(opts.filtType,'ht')
    
% initialize GPU variabes and wavelet filters
Vrec = zeros(m,n,opts.numMax,opts.cycleSpin,'single');
if gpuDeviceCount>0
    Vrec = gpuArray(Vrec);
    Volume = gpuArray(Volume);
end
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
    if gpuDeviceCount>0, Vavg = max(gather(mean(Vrec,4)),0);
    else, Vavg = max(mean(Vrec,4),0);
    end
        
elseif strcmp(opts.filtType,'median') % median filter version for salt and pepper noise
    % tau = 2*ones(1,3)*round(mean(tau(:))*15) + 1
    tau = [3 3 3];
    Vavg = (medfilt3(Volume,tau));
else
    error('unrecognized filtering');
end


% aggregate results, looping over each block
d2 = opts.matchSize;
U = zeros(m+d2,n+d2); % pad for block matches that wrap around
W = zeros(m+d2,n+d2);
out.Z = zeros(m,n);
D = getBM3DAggWeights(3,1);

mask = linspace(.1,.9,d);
mask = sin(pi*mask);
mask = mask'*mask;

% aggregate the blocks back into image and keep track of weights
for y = 1:M
    for x = 1:N
        Nm = numel(S{y,x});
        
        blocks = Vavg((y-1)*d+1:y*d,(x-1)*d+1:x*d,1:Nm);
        weight = 1/(D(blocks)+eps);

        for i = 1:Nm
            % get indices and shift according to reference location
            [a,b] = ind2sub([d2,d2],S{y,x}(i));            
            a = a + (y-1)*d;
            b = b + (x-1)*d;
            
            % aggregate into denoised image and add in weights
            U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight*blocks(:,:,i).*mask;
            W(a:a+d-1,b:b+d-1) = W(a:a+d-1,b:b+d-1) + weight.*mask;            
        end
        out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = Nm;
    end
end

U = U(d2/2+1:m+d2/2,d2/2+1:n+d2/2); % delete the padding
W = W(d2/2+1:m+d2/2,d2/2+1:n+d2/2);