N = 8;K = 16;
B = 512/8;


V = rand(N,N,K,B,B);

tic;
V2 = dct(V,[],1);
toc;

V1 = rand(N,N,K,B,B*2);
tic;
V3 = dct(V1,[],1);
toc;