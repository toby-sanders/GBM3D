d = 256;
SNR = 10;
levels = 3;
order = 2;
tau = .25;



wname = ['db',num2str(order)];
% P = phantom(d);
P = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
PO = P;
P = add_Wnoise(P,SNR);
tic;
C = myWavDec2D(wname,levels,P);
for i = 1:levels
    for j = 1:3
        C{i,j} = C{i,j}.*(abs(C{i,j})>tau);
    end
end

P2 = myWavRec2D(wname,levels,C);



shiftCnt = 3;
Pall = zeros(d,d,shiftCnt);
for ii = 0:shiftCnt-1
    C = myWavDec2D(wname,levels,circshift(P,[ii,ii]));
    for i = 1:levels
        for j = 1:3
            C{i,j} = C{i,j}.*(abs(C{i,j})>tau);
        end
    end

    Pall(:,:,ii+1) = circshift(myWavRec2D(wname,levels,C),[-ii,-ii]);
end
P3 = mean(Pall,3);
    




fprintf('time elapsed: %g\n',toc);
fprintf('noisy error: %g\n',myrel(P,PO));
fprintf('denoised error: %g\n',myrel(P2,PO));
fprintf('cyclespin error: %g\n',myrel(P3,PO));

figure(141);
subplot(2,2,1);imagesc(P,[0 1]);title('noisy');
subplot(2,2,2);imagesc(P2,[0 1]);title('denoised');
subplot(2,2,3);imagesc(P3,[0 1]);title('cycle spin');