clear;
d = 128;
SNR = 1e1;
levels = 3;
order = 1;
tau = 0.0;
rng(2021);

addpath('waveDecAndRec');
wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load(['waveletFilters',filesep,'db1Filters']);
elseif strcmp(wname,'db2'), load(['waveletFilters',filesep,'db2Filters']);
elseif strcmp(wname,'db3'), load(['waveletFilters',filesep,'db3Filters']);
end
loF = fft(Lo_d',d);
hiF = fft(Hi_d',d);

x = zeros(d,1);
% x = sin(1.5*pi*linspace(0,1,d)');
% x(d/4:3*d/4) = 1;
% x0 = x;

x(1:3*d/4) = 1;
% x(:) = 1;
x = sin(1.7*pi*linspace(0,1,d)');

% x = add_Wnoise(x,SNR);
tic;
C = myWavDec1D(wname,levels,x);
toc;
for j = 1:levels
    C{j} = C{j}.*(abs(C{j})>tau);
end
tic;
U2 = myWavRec1D(wname,levels,C);
toc;
U2 = U2(1:d);

x2 = zeros(2*d,1);
x2(1:d) = x;
x2(d+1:end) = linspace(x(end),x(1),d);



shiftCnt = 10;
    
Uall = zeros(d,shiftCnt);
for i = 1:shiftCnt
    C2 = myWavDec1D(wname,levels,circshift(x2,i-1));
    for j = 1:levels
        % C2{j}= max(abs(C2{j})-tau,0).*sign(C2{j});
        C2{j}= C2{j}.*(abs(C2{j})>tau);
    end
U3 = circshift(myWavRec1D(wname,levels,C2),-i+1);
Uall(:,i) = U3(1:d);
end
U3 = mean(Uall,2);

C = myWavDec1D_Mar_2021(wname,levels,x);
U4 = myWavRec1D(wname,levels,C);
U5 = myWavRec1DFast(wname,levels,C);
% 
% shiftCnt = 10;
% Uall = zeros(d,shiftCnt);
% for i = 1:shiftCnt
%     C = myWavDec(wname,levels,x);
%     for j = 1:levels
%         tmp = C{j}(1);
%         C{j}= max(abs(C{j})-tau,0).*sign(C{j});
%         C{j}(1) = tmp;
%     end
% 



myrel(x,U3)
myrel(x,U4)
myrel(x,U5)

for i = 1:levels+1
    subplot(3,3,i);plot(C{i});title(sprintf('level %i',i));
end
subplot(3,3,7);plot(U3);title('rec padded');
subplot(3,3,8);plot(U2);title('rec');
subplot(3,3,9);plot(x);title('true');

