function S = matchBlocksGPU(I,sigma,opts)

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
[m,n] = size(I);
tol = tol*d^2;

% current strategy is to tile the image with reference blocks
% a better strategy may involve keeping track to the pixels that get
% matched the most and searching for reference blocks that are missed most
M = floor(m/d);
N = floor(n/d);
S = cell(M,N);
I = gpuArray(single(I));

% precompute a basic 2D wavelet filtered image to improve the matching
if sigma>30/255
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
    I = Iavg;
end






% precompute norm of each dxd block in the image (blockNorms)
g = gpuArray(zeros(m,n,'single'));
g(1:d,1:d) = 1;
blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));

% loop over each reference block
FI = fft2(I);
parfor i = 1:M
    for j = 1:N
        % tic;
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        blockDistances = real(ifft2(conj(fft2(xR,m,n)).*(FI)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) + blockNorms ...
            - 2*blockDistances)/blockNorms((i-1)*d+1,(j-1)*d+1);
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
        S{i,j} = uint32(gather(S{i,j}));
        % toc;
%         tmp = zeros(m,n);
%         tmp(S{i,j}) = 1;
%         figure(123);imagesc(tmp);
%         pause;
    end
end
