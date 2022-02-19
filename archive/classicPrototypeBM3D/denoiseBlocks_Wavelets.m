function [U,out] = denoiseBlocks(I,S,d,sigma,opts)

% collaborative image denoising based on matched blocks
% inputs:
%   I - noisy grayscale image, values in [0,1]
%   S - cell matrix containing indices of matched blocks
%   d - block size
%   sigma - noise standard deviation
%   opts - options

% output: U - denoised image

tau = 3.5*sigma
if ~isfield(opts,'numMax')
    opts.numMax = nan;
end

% parameters for denoising
order = 3;
levels = 3;
wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

% build all 8 orthogonal 3D filters from the two 1D filters
Psi = zeros(N,N,N,8);
Psi(:,:,:,1) =  Lo_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,2) =  Hi_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,3) =  Lo_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,4) =  Lo_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,5) =  Hi_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,6) =  Hi_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,7) =  Lo_d'*Hi_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,8) =  Hi_d'*Hi_d.*reshape(Hi_d,1,1,N);

% initialize other vars
[A,B] = size(S);
[m,n] = size(I);
U = zeros(m,n);
out.K = zeros(m,n);
out.Z = zeros(m,n);
Ipad = zeros(2*m,2*n);
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
pp = waitbar(0,'denoising');

% loop over each reference block
for y = 1:A
    for x = 1:B
       % fprintf('denoising block %i of %i\n',(y-1)*B + x,A*B);
        waitbar(((y-1)*B + x)/A/B,pp,sprintf('denoising block %i of %i\n',(y-1)*B + x,A*B));
        Nm = numel(S{y,x});
        if Nm>1
            % gather all the matched blocks into 3D array
            blocks = zeros(d,d,Nm); 
            for i = 1:Nm
                [a,b] = ind2sub([m,n],S{y,x}(i));
                blocks(:,:,i) = Ipad(a:a+d-1,b:b+d-1);
            end

            if ceil(log(Nm)/log(2))~=log(Nm)/log(2)
                padSize = 2^(ceil(log(Nm)/log(2))) - Nm;
                blocks = cat(3,blocks,blocks(:,:,1:padSize));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % hard wavelet thresholding
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %  blocks = [blocks,fliplr(blocks);flipud(blocks),fliplr(flipud(blocks))];
            C = myWavDec3D_preallocated(Psi,levels,blocks);
            for i = 1:levels
                for j = 1:7
                    C{i,j} = C{i,j}.*(abs(C{i,j})>tau);
                end
            end
            blocks = myWavRec3D_preallocated(Psi,levels,C);
            
            % aggregate the blocks back into image and add weights
            for i = 1:Nm
                [a,b] = ind2sub([m,n],S{y,x}(i));
                if a+d<=m & b+d<=n
                    U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + blocks(:,:,i);
                    out.K(a:a+d-1,b:b+d-1) = out.K(a:a+d-1,b:b+d-1) + 1;
                end
            end
            
        end
        out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = Nm;
    end
end

% reweight image by number of times each pixel is visited
U = U./out.K;
close(pp);