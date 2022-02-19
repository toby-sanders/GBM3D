function [S,Volume] = matchBlocksCPU(I,sigma,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% match blocks using convolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ||x-y||_2^2 = ||x||^2 + ||y||^2 - 2<x,y>
% the first two norms can be precomputed instantly on every block using a
% cross-correlation
% the last term can be compute for each block x by cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the returned variable S is a cell variable where each cell contains the
% list of indices matched to that reference block. The index corresponds to
% the top left pixel in the image


% empirically determined threshold
if size(sigma,1)>1, sigma = mean(sqrt(sigma(:))); end
tol = .005*sigma^2;

% initialize constants
d = opts.blockSize;
d2 = opts.matchSize;
[m,n] = size(I);
tol = tol*d^2;
M = floor(m/d); % number of tiles in y-coordinate
N = floor(n/d); % number of tiles in x-coordinate
I = (single(I));

% precompute a basic 2D wavelet filtered image to improve the matching
% this is removed for now, will need a fix for building the volume with the
% original image 

tau = 0;% 2*sigma;
[Psi,Psi2] = getWaveFilters2Dbior(15);

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
    Iavg = Iavg + circshift(myWavRec2DGPUfast(Psi2,levels,C),[1-k,1-k])/shiftCnt;
end
Id = gather(Iavg);



% make padded image
Ipad = zeros(2*m,2*n,'single'); 
Ipad(1:m,1:n) = Id; Ipad(1:m,n+1:end) = Id; 
Ipad(m+1:end,1:n) = Id;Ipad(m+1:end,n+1:end) = Id;
Ipad = Ipad(1:m+d2,1:n+d2);
Ipad = circshift(Ipad,[d2/2,d2/2]);

Ipad2 = zeros(2*m,2*n,'single'); 
Ipad2(1:m,1:n) = I; Ipad2(1:m,n+1:end) = I; 
Ipad2(m+1:end,1:n) = I;Ipad2(m+1:end,n+1:end) = I;
Ipad2 = Ipad2(1:m+d2,1:n+d2);
Ipad2 = circshift(Ipad2,[d2/2,d2/2]);

% precompute norm of each dxd block in the image (blockNorms)
g = ones(d,'single');
blockNorms = real(ifft2(conj(fft2(g,m,n)).*(fft2(Id.^2))));
blockNormsPad = real(ifft2(conj(fft2(g,m+d2,n+d2)).*fft2(Ipad.^2)));

% initialize volume
Volume = zeros(m,n,opts.numMax,'single'); 
blocks = zeros(d,d,opts.numMax,'single');
S = cell(M,N);
% loop over each reference block
for i = 1:M
    indI = (i-1)*d;
    for j = 1:N
        indJ = (j-1)*d;
        
        xR = I((i-1)*d+1:i*d,(j-1)*d+1:j*d); % reference block
        Itmp = Ipad(indI+1:indI+d+d2,indJ+1:indJ+d+d2); % local image patch
        Itmp2 = Ipad2(indI+1:indI+d+d2,indJ+1:indJ+d+d2); % local image patch

        % cross-correlation of ref. block with local image patch
        blockDistances = real(ifft2(conj(fft2(xR,d+d2,d+d2)).*fft2(Itmp)));       
        blockDistances = (blockNorms((i-1)*d+1,(j-1)*d+1) +...
            blockNormsPad(indI+1:indI+d2,indJ+1:indJ+d2) ...
            - 2*blockDistances(1:d2,1:d2))/blockNorms((i-1)*d+1,(j-1)*d+1);
        
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
        % with real data due to numerical error
        S{i,j}(1) = d2^2/2 + d2/2 + 1;
        
        % add matched blocks to volume
        
        for k = 1:numel(S{i,j})
            [a,b] = ind2sub([d2,d2],S{i,j}(k));
            blocks(:,:,k) = Itmp2(a:a+d-1,b:b+d-1);
        end
        Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,1:numel(S{i,j})) = blocks;
        
        % pad the rest of the block in the volume with the reference block
        Volume((i-1)*d+1:i*d,(j-1)*d+1:j*d,k+1:opts.numMax) = repmat(xR,1,1,opts.numMax-k);
    end
end
