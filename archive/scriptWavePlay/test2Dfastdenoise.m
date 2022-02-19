d = 256;
SNR = 5;
order = 4;
levels = 4;
wname = 'sym';
testV = linspace(1,10,100);
myDir = pwd;

I = im2double(imread('cameraman.tif'));
[b,sigma] = add_Wnoise(I,SNR);
cd ../
Psi = getWaveFilters2D(wname,order);
cd(myDir);


recs = zeros(d,d,numel(testV));
for k = 1:numel(testV)
    C = myWavDec2DFFT(Psi,levels,b);

    tau = sigma*testV(k);
    ss = 0;ss2 = 0;
    for i = 1:levels
        for j = 1:3
            C{i,j} = C{i,j}.*(abs(C{i,j})>tau);        
        end
    end
    rec = gather(myWavRec2DGPUfast(Psi,levels,C));
    recs(:,:,k) = rec;
    ee(k) = myrel(rec,I);
end


myrel(rec,b)
myrel(rec,I)
[~,mm] = min(ee);
figure(114);colormap(gray);
subplot(2,2,1);imagesc(recs(:,:,mm),[0 1]);title('denoised');
subplot(2,2,2);imagesc(I,[0 1]);title('original');
subplot(2,2,3);imagesc(b,[0 1]);title('noisy');
subplot(2,2,4);plot(testV,ee);