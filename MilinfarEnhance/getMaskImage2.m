function U = getMaskImage2(I)



[m,n] = size(I);

d = 3;
d2 = 12;
step = 1;

% make padded image
Ipad = zeros(2*m,2*n,'single'); 
Ipad(1:m,1:n) = I; Ipad(1:m,n+1:end) = I; 
Ipad(m+1:end,1:n) = I;Ipad(m+1:end,n+1:end) = I;
Ipad = Ipad(1:m+d2+d,1:n+d2+d);
Ipad = circshift(Ipad,[d2/2+(d-1)/2,d2/2+(d-1)/2]);


% precompute norm of each dxd block in the image (blockNorms)
g = ones(d,'single');
blockNorms = real(ifft2(conj(fft2(g,m,n)).*(fft2(I.^2))));
blockNormsPad = real(ifft2(conj(fft2(g,size(Ipad,1),size(Ipad,2))).*fft2(Ipad.^2)));
U = zeros(m+d,n+d);
Nmatrix = zeros(m+d,n+d);
for i = 1:step:m
    for j = 1:step:n
        indI = d2/2+i-1;
        indJ = d2/2+j-1;
        
        xR = Ipad(indI+1:indI+d,indJ+1:indJ+d);
        Itmp = Ipad(i:i+d+d2-1,j:j+d+d2-1);
        
        % cross-correlation of ref. block with local image patch
        blockDistances = real(ifft2(conj(fft2(xR,d+d2,d+d2)).*fft2(Itmp))); 
        blockDistances = blockNormsPad(indI+1,indJ+1) +...
            blockNormsPad(i:i+d2-1,j:j+d2-1) ...
            - 2*blockDistances(1:d2,1:d2);
        V = mink(blockDistances(:),16);
        % U(i:i+d-1,j:j+d-1) = U(i:i+d-1,j:j+d-1)+ sum(blockDistances(:));
        % Nmatrix(i:i+d-1,j:j+d-1)  = Nmatrix(i:i+d-1,j:j+d-1)  + 1;
        U(i + (d-1)/2, j +(d-1)/2) = U(i + (d-1)/2, j +(d-1)/2)+ sum(V);% sum(blockDistances(:));
        Nmatrix(i + (d-1)/2, j +(d-1)/2)  = Nmatrix(i + (d-1)/2, j +(d-1)/2)  + 1;

    end
end
U = U./Nmatrix;
U = U((d+1)/2:(d-1)/2+m,(d+1)/2:(d-1)/2 + n);
