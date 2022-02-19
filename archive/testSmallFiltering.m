clear;
d = 64;

z = rand(d);
y = rand(d);


for i = 1:100
    tic;
    t = imfilter(z,y);
    toc;
end

gz = gpuArray(double(z));
gy = gpuArray(double(y));
for i =1:100
    tic;
    t = imfilter(gz,gy);
    toc;
end

gz = gpuArray(double(z));
gy = gpuArray(double(y));



z = rand(d,'single');
y = rand(d,'single');
for i = 1:100
    tic;
    t = ifft2(fft2(z).*conj(fft2(y)));
    toc;
end

gz = gpuArray((z));
gy = gpuArray((y));
for i = 1:100
    tic;
    t = ifft2(fft2(z).*conj(fft2(y)));
    toc;
end



