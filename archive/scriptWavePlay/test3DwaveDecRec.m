clear;
d = 128;
levels = 3;

[Psi,Psi2] = getWaveFilters3Dbior('db',15,'1');


U = rand(d,d,d);

C = myWavDec3DFFT(Psi,levels,U);

U2 = myWavRec3DGPUfast(Psi2,levels,C);

myrel(U,U2)

%%
norm(U(:))^2
s = 0;
for i = 1:numel(C)
    s = s + norm(C{i}(:))^2;
end
s