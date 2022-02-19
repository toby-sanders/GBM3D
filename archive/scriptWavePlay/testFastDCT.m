clear;
d = 8;
sZ = 16;
M = 64;
N = 64;



x = zeros(d,d,sZ,M,N);
x0 = phantom(d);
for i = 1:sZ
    for j = 1:M
        for k = 1:N
            x(:,:,i,j,k) = x0;
        end
    end
end

T = dct(eye(d),[],1);
T2 = dct(eye(d),[],2);
gX = gpuArray(x);
gT = gpuArray(T);

tic;
Tx = dct(x,[],1);
Tx = dct(Tx,[],2);
toc;

tic;
Tx2 = pagemtimes(T,x);
 Tx2 = pagemtimes(Tx2,T');
% Tx2 = pagetranspose(pagemtimes(T,'none',Tx2,'transpose'));
toc;

f = @()pagemtimes(gX,gT');
% @()pagetranspose(pagemtimes(gT,'none',gX,'transpose'));
t = gputimeit(f)

tic;
Tx3 = pagemtimes(gT,gX);
toc;

myrel(Tx,Tx2)