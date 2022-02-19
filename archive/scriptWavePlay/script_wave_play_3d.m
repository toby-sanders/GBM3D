d = 256;
SNR = 10;
levels = 3;
order = 3;
tau = .25;

wname = ['db',num2str(order)];
P = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
% P = phantom(d);
P = repmat(P,1,1,d/2);
P = add_Wnoise(P,SNR);

tic;
C = myWavDec3D(wname,levels,P);
for i = 1:levels
    for j = 1:7
        C{i,j} = C{i,j}.*(abs(C{i,j})>tau);
    end
end
P2 = myWavRec3D(wname,levels,C);
toc;


%%
figure(212);colormap(gray);
subplot(2,2,1);imagesc(P(:,:,1),[0 1]);
subplot(2,2,2);imagesc(P2(:,:,1),[0 1]);