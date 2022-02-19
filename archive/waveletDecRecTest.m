clear;
d = 128;
SNR = 5e17;
levels = 3;
order = 3;
tau = 0;
rng(2021);

wname = ['db',num2str(order)];

x = zeros(d,1);
x = sin(1.5*pi*linspace(0,1,d)');
x(d/4:3*d/4) = 1;
x0 = x;


c = myWavDec2(wname,levels,x);
rec = myWavRec2(wname,levels,c);

A1 = @(U)myWavDec2(wname,levels,U);
A2 = @(C)myWavRec2(wname,levels,C);

[flg,~,~] = check_D_Dt(A1,A2,[d,1,1]);
flg
norm(x)
norm(c)
norm(rec)