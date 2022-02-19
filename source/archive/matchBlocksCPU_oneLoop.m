function S = matchBlocksCPU(I,sigma,opts)

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
S = cell(M*N,1);
I = (single(I));

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
    I = gather(Iavg);
end






% precompute norm of each dxd block in the image (blockNorms)
g = (zeros(m,n,'single'));
g(1:d,1:d) = 1;

Ipad = zeros(2*m,2*n,'single'); % make padded image
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
Ipad = Ipad(1:m+d,1:n+d);
Ipad = circshift(Ipad,[d/2,d/2]);
blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));
blockNormsPad = real(circshift(ifft2(fft2(g,m+d,n+d).*fft2(Ipad.^2)),[-d+1,-d+1]));
% loop over each reference block
for ii = 1:M*N
    [i,j] = ind2sub([M,N],ii);
    % tic;
    xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
    tmp = Ipad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
    blockDistances = real(ifft2(conj(fft2(xR,2*d,2*d)).*fft2(tmp)));       
    blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) +...
        blockNormsPad((i-1)*d+1:(i)*d,(j-1)*d+1:(j)*d) ...
        - 2*blockDistances(1:d,1:d))/blockNorms((i-1)*d+1,(j-1)*d+1);

    if opts.numMax == opts.numMin
        [~,S{ii}] = mink(blockDistances(:),opts.numMax);
    else
        S{ii} = find(blockDistances<tol);

        % don't allow for too few or too many blocks in the matching
        if numel(S{ii})>opts.numMax
            [~,S{ii}] = mink(blockDistances(:),opts.numMax);
        elseif numel(S{ii})<opts.numMin
            [~,S{ii}] = mink(blockDistances(:),opts.numMin);
            % S{i,j} = [S{i,j};repmat(S{i,j}(2),opts.numMax-opts.numMin,1)];
        end       
    end
    S{ii} = uint32((S{ii}));
end
S = reshape(S,M,N);
