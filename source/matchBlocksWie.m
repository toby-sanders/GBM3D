function [S,V1,V2] = matchBlocksWie(I,In,sigma,opts)

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
if size(sigma,1)>1, sigma = mean(sqrt(sigma(:))); end
tol = .005*sigma^2;

d = opts.blockSizeWie;
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

Ipad = zeros(2*m,2*n,'single'); % make padded image
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
Ipad = Ipad(1:m+d2,1:n+d2);
Ipad = circshift(Ipad,[d2/2,d2/2]);

Ipad2 = zeros(2*m,2*n,'single'); % make padded image
Ipad2(1:m,1:n) = In; Ipad2(1:m,n+1:end) = In; 
Ipad2(m+1:end,1:n) = In;Ipad2(m+1:end,n+1:end) = In;
Ipad2 = Ipad2(1:m+d2,1:n+d2);
Ipad2 = circshift(Ipad2,[d2/2,d2/2]);

% blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));
% blockNormsPad = real(circshift(ifft2(fft2(g,m+d,n+d).*fft2(Ipad.^2)),[-d+1,-d+1]));
blockNorms = real(ifft2(conj(fft2(g)).*(fft2(I.^2))));
blockNormsPad = real(ifft2(conj(fft2(g,m+d2,n+d2)).*fft2(Ipad.^2)));
V1 = zeros(d,d,opts.numMax,M,N,'single'); % initialize volume
V2 = V1;

% loop over each reference block
for i = 1:M
    indI = (i-1)*d;
    for j = 1:N
        indJ = (j-1)*d;
        
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        Itmp = Ipad(indI+1:indI+d+d2,indJ+1:indJ+d+d2);
        Itmp2 = Ipad2(indI+1:indI+d+d2,indJ+1:indJ+d+d2);
        % Itmp = Ipad((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
        % Itmp2 = Ipad2((i-1)*d+1:(i+1)*d,(j-1)*d+1:(j+1)*d);
        blockDistances = real(ifft2(conj(fft2(xR,d+d2,d+d2)).*fft2(Itmp)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) +...
            blockNormsPad(indI+1:indI+d2,indJ+1:indJ+d2) ...
            - 2*blockDistances(1:d2,1:d2))/blockNorms((i-1)*d+1,(j-1)*d+1);
        
        [~,S{i,j}] = mink(blockDistances(:),opts.numMax);
        
        S{i,j} = uint32((S{i,j}));
        % ensure the reference block is first. This can be off sometimes
        % with real data due to numerical error
        S{i,j}(1) = d2^2/2 + d2/2 + 1;
        
        % add matched blocks to volume
        for k = 1:numel(S{i,j})
            [a,b] = ind2sub([d2,d2],S{i,j}(k));
            V1(:,:,k,i,j) = Itmp(a:a+d-1,b:b+d-1);
            V2(:,:,k,i,j) = Itmp2(a:a+d-1,b:b+d-1);
        end

    end
end
