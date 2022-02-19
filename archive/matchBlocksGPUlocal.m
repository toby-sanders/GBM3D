function S = matchBlocksGPU(I,sigma,opts)

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

% hardcoded constants to test
if ~isfield(opts,'blockSize')
    d = 8;
else
    d = opts.blockSize;
end
if ~isfield(opts,'numMin')
    opts.numMin = 32;
end
if ~isfield(opts,'numMax')
    opts.numMax = max(256,opts.numMin);
end
% C = 2.2605;
C = .005;
tol = C*sigma^2;


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

% Ipad = zeros(m+2*d,n+2*d);
% Ipad(d+1:m+d , d+1:n+d) = I;
Ipad = zeros(2*m,2*n,'gpuArray'); % make padded image
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
Ipad = Ipad(1:m+d,1:n+d);
Ipad = circshift(Ipad,[d/2,d/2]);
blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));
blockNormsPad = real(circshift(ifft2(fft2(g,m+d,n+d).*fft2(Ipad.^2)),[-d+1,-d+1]));
% loop over each reference block
for i = 1:M
    for j = 1:N
        % tic;
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        tmp = Ipad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
        blockDistances0 = real(ifft2(conj(fft2(xR,2*d,2*d)).*fft2(tmp)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) +...
            blockNormsPad((i-1)*d+1:(i)*d,(j-1)*d+1:(j)*d) ...
            - 2*blockDistances0(1:d,1:d))/blockNorms((i-1)*d+1,(j-1)*d+1);
        % blockDistances = blockDistances(1:2*d,1:2*d);
        S{i,j} = find(blockDistances<tol);
%         t2 = blockNormsPad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
%         figure(12);
%         subplot(2,2,1);imagesc(blockDistances0);colorbar;
%         subplot(2,2,2);imagesc(tmp);colorbar;
%         subplot(2,2,3);imagesc(blockDistances);colorbar;
%         subplot(2,2,4);imagesc(t2);colorbar;
        
        % don't allow for too few or too many blocks in the matching
        if numel(S{i,j})>opts.numMax
            [~,S{i,j}] = mink(blockDistances(:),opts.numMax);
        elseif numel(S{i,j})<opts.numMin
            [~,S{i,j}] = mink(blockDistances(:),opts.numMin);
            % S{i,j} = [S{i,j};repmat(S{i,j}(2),opts.numMax-opts.numMin,1)];
        end             
        S{i,j} = uint32(gather(S{i,j}));
        % toc;
%         tmp = zeros(m,n);
%         tmp(S{i,j}) = 1;
%         figure(123);imagesc(tmp);
%         pause;
    end
end
