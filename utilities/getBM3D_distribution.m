function [R,d0,k] = getBM3D_distribution(I,matchSpins,k)

% set up the parameters for distributing the workload for the GBM3D
% processing to run in parallel

% Written by Toby Sanders @Lickenbrock Tech
% 3-17-2021


d0 = 256+k;

[m,n] = size(I);
Xcnt = ceil((n-k)/(d0-k));
Ycnt = ceil((m-k)/(d0-k));
Tcnt = Xcnt*Ycnt*numel(matchSpins);
R = zeros(Tcnt,5);

cntI = 0;
cntJ = 0;
cnt = 0;

ramp = repmat(linspace(.01,1,k),256,1);

for i = 1:Ycnt
    indY = (i-1)*(d0-k);
    for j = 1:Xcnt
        indX = (j-1)*(d0-k);
        for mS = 1:numel(matchSpins)
            cnt = cnt+1;
            R(cnt,1) = indY+1;
            R(cnt,2) = indY+d0;
            R(cnt,3) = indX+1;
            R(cnt,4) = indX+d0;
            R(cnt,5) = matchSpins(mS);
            if R(cnt,4)>=n
               R(cnt,3) = min(n-ceil((n-R(cnt,3)+1)/k)*k+1,n-2*k+1);
            end
            if R(cnt,2)>=m
               R(cnt,1) = min(m-ceil((m-R(cnt,1)+1)/k)*k+1,m-2*k+1);
            end
        end
        cntI = cntI+d0-32;
        cntJ = cntJ+d0-32;
    end
end

R(:,2) = min(R(:,2),m);
R(:,4) = min(R(:,4),n);