function [U,out] = denoiseBlocks(I,S,d,sigma,opts)

% collaborative image denoising based on matched blocks
% inputs:
%   I - noisy grayscale image, values in [0,1]
%   S - cell matrix containing indices of matched blocks
%   d - block size
%   sigma - noise standard deviation
%   opts - options

% output: U - denoised image

if ~isfield(opts,'tau')
    tau = 3*sigma;
else
    tau = opts.tau;
end
if ~isfield(opts,'numMax')
    opts.numMax = nan;
end

% parameters for HOTV denoising code
hopts.order = 2;
hopts.levels = 1;
hopts.iter = 75;
hopts.wrap_shrink = false;

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
            
            % denoise the 3D array with hard thresholding
            if Nm == opts.numMax, hopts.wname = 'db1';
            else, hopts.wname = 'db2';
            end
            hopts.levels = 4;
            rec = my3DwavDenoise(blocks,tau,hopts);
            
            % aggregate the blocks back into image and add weights
            for i = 1:Nm
                [a,b] = ind2sub([m,n],S{y,x}(i));
                if a+d<=m & b+d<=n
                    U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + rec(:,:,i);
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