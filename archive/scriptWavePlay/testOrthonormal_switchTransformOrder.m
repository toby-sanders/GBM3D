clear;
Psi = single(reshape([-sqrt(2)/2,sqrt(2)/2],1,1,2));
Phi = abs(Psi);



V = ones(8,8,16);
V1 = dct(V,[],1);
V1 = dct(V1,[],2);
V1 = myWavDec1Don3D(Psi,Phi,2,V1);

V2 = myWavRec1Don3D(Psi,Phi,2,V1);
V2 = idct(V2,[],2);
V2 = idct(V2,[],1);

V3 = idct(V1,[],2);
V3 = idct(V3,[],1);
V3 = myWavRec1Don3D(Psi,Phi,2,V3);

V4 = myWavRec1Don3D(Psi,Phi,2,V1);
V4 = idct(V4,[],1);
V4 = idct(V4,[],2);


myrel(V2,V)
myrel(V3,V)
myrel(V4,V)