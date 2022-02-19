function [U,out] = denoiseAndAggregateGPU(Volume,I,S,tau,Psi,Psi2,opts)

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


% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
[m,n] = size(I);
d = opts.blockSize;


% aggregate results, looping over each block
U = zeros(m+d,n+d); % pad for block matches that wrap around
out.K = zeros(m+d,n+d);
out.Z = zeros(m,n);
D = getBM3DAggWeights(3,1);

Psi = gather(Psi);Psi2 = gather(Psi2);
W = linspace(.1,.9,d);
W = sin(pi*W);
W = W'*W;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard wavelet thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% opts.cycleSpin: number of spins in the wavelet cycle spinning
if strcmp(opts.filtType,'ht')
    for y = 1:M
        for x = 1:N
            xR = Volume((y-1)*d+1:y*d,(x-1)*d+1:x*d,:);
            % wavelet decomposition
            C = myWavDec3DFFT(Psi,opts.levels,xR);   

            % hard thresholding
            for i = 1:opts.levels 
                for j = 1:7
                    C{i,j} = C{i,j}.*(abs(C{i,j})>tau(i,j));
                end
            end

            % wavelet recontruction
            blocks = myWavRec3DGPUfast(Psi2,opts.levels,C);

             weight = 1/(D(blocks)+eps);
             Nm = numel(S{y,x});

            for i = 1:Nm
                [a,b] = ind2sub([d,d],S{y,x}(i));
                a = a + (y-3/2)*d;
                b = b + (x-3/2)*d;
                if a>0 & b>0
                U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight*blocks(:,:,i).*W;
                out.K(a:a+d-1,b:b+d-1) = out.K(a:a+d-1,b:b+d-1) + weight.*W;
                end
            end
            out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = Nm;
        end
    end
end


% reweight image by number of times each pixel is visited
U = U./out.K;
U = U(1:m,1:n); % delete the padding
out.K = out.K(1:m,1:n);