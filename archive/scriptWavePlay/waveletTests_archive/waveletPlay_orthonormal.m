clear;
d = 1280;
SNR = 10;
levels = 6;
order = 1;
tau = 0.3;
rng(2021);

wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
loF = fft(Lo_d',d);
hiF = fft(Hi_d',d);

x = zeros(d,1);
x = sin(1.5*pi*linspace(0,1,d)');
x(d/4:3*d/4) = 1;
x0 = x;
x = add_Wnoise(x,SNR);


cnt = 0;
for j = 1:levels
    cnt = cnt+d/2^j;
end
    

% U = my_wavelet_denoise_1D_ortho2(wname,levels,x,tau);
shiftCnt = 20;
Uall = zeros(d,shiftCnt);
for i = 1:shiftCnt
    C = myWavDec(wname,levels,circshift(x,i));
    for j = 1:levels
        % ind = mod(i,numel(C{j}))
        % tmp = C{j}(ind);
        C{j}= max(abs(C{j})-tau,0).*sign(C{j});
       %  C{j}(ind) = tmp;
    end
    Uall(:,i) = circshift(myWavRec(wname,levels,C),-i);
    
end
C = myWavDec(wname,levels,x);
for j = 1:levels
    tmp = C{j}(1);
    C{j}= max(abs(C{j})-tau,0).*sign(C{j});
    C{j}(1) = tmp;
end
U2 = myWavRec(wname,levels,C);
% tau = .4;
% mm = 16;
% [c,l] = wavedec(x,levels,wname);
% c(mm+1:end) = c(mm+1:end).*(abs(c(mm+1:end))>tau);
% wrec = waverec(c,l,wname);


hopts.mu = 3;
hopts.order = order;
hopts.mode = 'deconv';
h = zeros(d,1);
h(1) = 1;
rec = HOTV3D(h,x,[d,1,1],hopts);

figure(121);hold off;
subplot(2,2,1);hold off;
plot(x,'linewidth',2);title('signal');hold on;
plot(x0,'--','linewidth',1.5);hold off;
subplot(2,2,2);hold off;
plot(real(U2),'-','linewidth',2);title('recon ortho wavelets');hold on;
plot(x0,'--','linewidth',1.5);hold off;
subplot(2,2,3);hold off;
plot(rec,'linewidth',2);title('TV denoised');hold on;
plot(x0,'--','linewidth',1.5);hold off;
subplot(2,2,4);hold off;
plot(mean(Uall,2),'linewidth',2);title('cycle spinning');hold on;
plot(x0,'--','linewidth',1.5);hold off;
hold off;

myrel(U2,x0)
myrel(mean(Uall,2),x0)
myrel(rec,x0)
myrel(x,x0)