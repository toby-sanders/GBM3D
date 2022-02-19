function [R,d0] = getBM3D_distribution(I,matchSpins)

d0 = 256;
[m,n] = size(I);

numPatches = m*n/d0^2;
numCases = numPatches*numel(matchSpins);

R = zeros(numCases,3);
cnt = 0;
for i = 1:m/d0
    for j = 1:n/d0
        for k = 1:numel(matchSpins)
            cnt = cnt+1;
            R(cnt,1) = i;
            R(cnt,2) = j;
            R(cnt,3) = matchSpins(k);
        end
    end
end
