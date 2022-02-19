clear;
d = 256;
SNR = 5e0;
order = 15;
levels = 3;
wname = 'bior';
testV = 10;% linspace(10,100,100);
myDir = pwd;

I = phantom(d);
I = I(:,d/2);
% I = im2double(imread('cameraman.tif'));
[b,sigma] = add_Wnoise(I,SNR);
cd ../waveletFilters
load([wname,num2str(order),'Filters']);
cd(myDir);

if exist('Hi_d')
    HiD = Hi_d;
    LoD = Lo_d;
    HiR = Hi_d;
    LoR = Lo_d;
else
    LoR = fliplr(LoR);
    HiR = fliplr(HiR);
end



C = myWavDec1DBior(LoD,HiD,levels,b);

tau = sigma*testV;
ss = 0;ss2 = 0;
nn = 0;
for i = 1:levels
    nn = nn + norm(C{i})^2;
    % for j = 1:3
        C{i} = C{i}.*(abs(C{i})>tau);        
    % end
end
nn = nn + norm(C{end})^2;
rec = myWavRec1DBior(LoR,HiR,levels,C);


nn
norm(I)^2
norm(rec(:))^2
figure(114);colormap(gray);
subplot(2,2,1);plot(rec);title('denoised');
subplot(2,2,2);plot(I);title('original');
subplot(2,2,3);plot(b);title('noisy');
% subplot(2,2,4);plot(testV,ee);