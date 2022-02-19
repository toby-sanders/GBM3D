function [S] = matchBlocksGPU(I,sigma,opts)

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
tau = 2.2605;
tol = tau*sigma^2;
d = opts.blockSize;
d2 = opts.matchSize;

[m,n] = size(I);
tol = tol*d^2;

% current strategy is to tile the image with reference blocks
% a better strategy may involve keeping track to the pixels that get
% matched the most and searching for reference blocks that are missed most
M = floor(m/d);
N = floor(n/d);
S = cell(M,N);
I = (single(I));

% precompute norm of each dxd block in the image (blockNorms)
g = (zeros(m,n,'single'));
g(1:d,1:d) = 1;
blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));

% loop over each reference block
FI = fft2(I);

for i = 1:M
    for j = 1:N
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        blockDistances = real(ifft2(conj(fft2(xR,m,n)).*(FI)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) + blockNorms - 2*blockDistances);
        
        S{i,j} = find(blockDistances<tol);
        
        % don't allow for too few or too many blocks in the matching
        if numel(S{i,j})>opts.numMax
            [~,S{i,j}] = mink(blockDistances(:),opts.numMax);
        elseif numel(S{i,j})<opts.numMin
            [~,S{i,j}] = mink(blockDistances(:),opts.numMin);
        end             
        S{i,j} = uint32((S{i,j}));
        
    end
end
