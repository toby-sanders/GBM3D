function [S,Volume] = matchBlocksV(I,sigma,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% match blocks using convolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ||x-y||_2^2 = ||x||^2 + ||y||^2 - 2<x,y>
% the first two norms can be precomputed instantly on every block using a
% convolution
% the last term can be compute for each block x by convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the returned variable S is a cell variable where each cell contains the
% list of indices matched to that reference block. The index corresponds to
% the top left pixel in the image


% empirically determined threshold
tol = .005*sigma^2;

d = opts.blockSize;
[m,n,K] = size(I);
if mod(K,2)~=1, error('input volume should have odd number of frames'); end
tol = tol*d^2;
ref = (K+1)/2;

% current strategy is to tile the image with reference blocks
% a better strategy may involve keeping track to the pixels that get
% matched the most and searching for reference blocks that are missed most
M = floor(m/d);
N = floor(n/d);
S = cell(M,N);
I = (single(I));

% precompute a basic 2D wavelet filtered image to improve the matching
% this is removed for now, will need a fix for building the volume with the
% original image 
if sigma>1e20% 30/255
    tau = 2*sigma;
    Psi = getWaveFilters2D('sym',4);
    Iavg = gpuArray(zeros(size(I),'single'));
    shiftCnt = 5;
    levels = 3;
    for k = 1:shiftCnt
        C = myWavDec2DFFT(Psi,levels,circshift(I,[k-1,k-1]));
        for i = 1:levels
            for j = 1:3
                C{i,j} = C{i,j}.*(abs(C{i,j})>tau);        
            end
        end
        Iavg = Iavg + circshift(myWavRec2DGPUfast(Psi,levels,C),[1-k,1-k])/shiftCnt;
    end
    I = gather(Iavg);
end






% precompute norm of each dxd block in the image (blockNorms)
g = (zeros(m,n,'single'));
g(1:d,1:d) = 1;

Ipad = zeros(2*m,2*n,K,'single'); % make padded image
Ipad(1:m,1:n,:) = I; Ipad(1:m,n+1:end,:) = I; 
Ipad(m+1:end,1:n,:) = I;Ipad(m+1:end,n+1:end,:) = I;
Ipad = Ipad(1:m+d,1:n+d,:);
Ipad = circshift(Ipad,[d/2,d/2,0]);
% blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));
% blockNormsPad = real(circshift(ifft2(fft2(g,m+d,n+d).*fft2(Ipad.^2)),[-d+1,-d+1]));
blockNorms = real(ifft2(conj(fft2(g)).*(fft2(I.^2))));
blockNormsPad = real(ifft2(conj(fft2(g,m+d,n+d)).*fft2(Ipad.^2)));
Volume = zeros(m,n,opts.numMax,'single'); % initialize volume
blocks = zeros(d,d,opts.numMax,'single');
% loop over each reference block
for i = 1:M
    for j = 1:N
        % tic;
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d,ref); % reference block
        Itmp = Ipad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d,:);
        blockDistances = real(ifft2(conj(fft2(xR,2*d,2*d)).*fft2(Itmp)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1,ref) +...
            blockNormsPad((i-1)*d+1:(i)*d,(j-1)*d+1:(j)*d,:) ...
            - 2*blockDistances(1:d,1:d,:))/blockNorms((i-1)*d+1,(j-1)*d+1,ref);
        
        if opts.numMax == opts.numMin
            [~,S{i,j}] = mink(blockDistances(:),opts.numMax);
        else
            S{i,j} = find(blockDistances<tol);

            % don't allow for too few or too many blocks in the matching
            if numel(S{i,j})>opts.numMax
                [~,S{i,j}] = mink(blockDistances(:),opts.numMax);
            elseif numel(S{i,j})<opts.numMin
                [~,S{i,j}] = mink(blockDistances(:),opts.numMin);
                % S{i,j} = [S{i,j};repmat(S{i,j}(2),opts.numMax-opts.numMin,1)];
            end       
        end
        S{i,j} = uint32((S{i,j}));
        % ensure the reference block is first. This can be off sometimes
        % due to numerical error
        S{i,j}(1) = d^2*K/2 + d/2 + 1;% d*d*(K-1)/2 + d*d/2 + d/2+1; 
        
        % add matched blocks to volume
        for k = 1:numel(S{i,j})
            [a,b,c] = ind2sub([d,d,K],S{i,j}(k));
            blocks(:,:,k) = Itmp(a:a+d-1,b:b+d-1,c);
        end
        Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,1:numel(S{i,j})) = blocks;
        
        % pad the rest of the block in the volume with the reference block
        Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,k+1:opts.numMax) = repmat(xR,1,1,opts.numMax-k);
    end
end
