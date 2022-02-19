clear;
d = 32;
sZ = 16;
SNR = 1000;
profile = 'fast';
levels = 2;

I = phantom(d);
[bb,sigma] = add_Wnoise(I,SNR);


opts = setBM3Dopts(profile);
d2 = opts.matchSize;

[S,V] = matchBlocksCPU(bb,sigma,opts);
P = zeros(numel(V),d^2);
cnt = 0;
for i = 1:numel(S)
    for i2 = 1:size(S{i})
        [a,b] = ind2sub([d2,d2],S{i}(i2));
        a = a - d2/2;
        b = b - d2/2;
        for k = 1:sZ
            for j = 1:sZ
                cnt = cnt+1;
                a2 = mod(a+(j-1)-1,d2) +1;
                b2 = mod(b+(k-1)-1,d2) +1;
                
                tmp = sub2ind([d2,d2],a2,b2);
                P(cnt,tmp) = 1;
                
            end
        end
    end
end

V2 = P*bb(:);
V3 = V;
cnt = 0;
for xx = 1:2
for yy = 1:2
for i1 = 1:16
for i2 = 1:sZ
for i3 = 1:sZ
    cnt = cnt+1;
    indz = i1;
    indx = (xx-1)*sZ + i2;
    indy = (yy-1)*sZ + i3;
    V3(indy,indx,indz) = V2(cnt);
    
end
end
end
end
end
[Psi,Psi2] = getWaveFilters3Dbior;
E = zeros(d,d,16);
T = zeros(d*d*16);
cnt1 = 0;

for xx = 1:2
for yy = 1:2
for i1 = 1:16
for i2 = 1:sZ
for i3 = 1:sZ
    indz = i1;
    indx = (xx-1)*sZ + i2;
    indy = (yy-1)*sZ + i3;
    E(indy,indx,indz) = 1;
    C = myWavDec3DFFT(Psi,levels,E);
    E(indy,indx,indz) = 0;
    
    C2 = E;
    cnt = 0;
    for j = 1:levels
        for k = 1:7
        C2(cnt+1:cnt+numel(C{j,k})) = C{j,k}(:);
        cnt = cnt + numel(C{j,k});
        end
    end
    C2(cnt+1:end) = C{levels+1,1};
    cnt1 = cnt1+1;
    T(cnt1,:) = reshape(C2,1,d^2*16);
end
end
end
end
end
