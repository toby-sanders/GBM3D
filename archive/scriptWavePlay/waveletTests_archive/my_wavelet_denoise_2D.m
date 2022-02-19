function U = my_wavelet_trans_2D(wname,levels,U)

% Written by Toby Sanders

% these are designed for 1-D problems
[p,q,r] = size(U);
% [Lo_d,Hi_d] = wfilters(wname);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end

L = fft(Lo_d',p);
H = fft(Hi_d',p);

m = numel(Lo_d);

C = zeros(p,levels+1);
FU = fft(U);
for i = 1:levels
    C(:,i) = ifft(FU.*H);
    FU = FU.*L;
end
C(:,levels+1) = ifft(FU);

tau = .1;
C(:,1:levels) = max(abs(C(:,1:levels))-tau,0).*sign(C(:,1:levels));



figure(223);hold off;
for i = 1:levels+1
    figure(223);
    plot(C(:,i));hold on;
end


C = fft(C);
U(:) = 0;
LHF = H;
for i = 1:levels
    U = U + ifft(C(:,i).*conj(LHF))/(2^i);
    LHF = LHF.*L;
end
U = U + ifft(C(:,levels+1).*conj(L.^levels))/(2^levels);



