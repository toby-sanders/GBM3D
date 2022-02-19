function U = myWavRec3D(wname,levels,C)

% 3-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
[p,q,r] = size(C{1});
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);
U = zeros(2*p,2*q,2*r); % store restored signal into U

% build all 8 orthogonal 3D filters from the two 1D filters
Psi = zeros(N,N,N,8);
Psi(:,:,:,1) =  Lo_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,2) =  Hi_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,3) =  Lo_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,4) =  Lo_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,5) =  Hi_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,6) =  Hi_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,7) =  Lo_d'*Hi_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,8) =  Hi_d'*Hi_d.*reshape(Hi_d,1,1,N);


Cup = cell(levels+1,7);
% HPF the upsampled detail coefficients
for i = 1:levels
    for j = 1:7
        Cup{i,j} = zeros(2*size(C{i,j},1),2*size(C{i,j},2),2*size(C{i,j},3));
        Cup{i,j}(1:2:end,1:2:end,1:2:end) = C{i,j};
        Cup{i,j} = circshift(imfilter(Cup{i,j},Psi(:,:,:,j+1),'corr','circular'),[-N/2+1,-N/2+1,-N/2+1]);
    end
end
Cup{levels+1,1} = C{levels+1,1};

% iteratively LPF and upsample
for i = 1:levels
    for j = 1:7
        U = U + Cup{i,j}; % next level coefficients added
    end
    % upsample coefficients and LPF
    for k = i+1:levels+1
        tmp = zeros(2*size(Cup{k,j},1),2*size(Cup{k,j},2),2*size(Cup{k,j},3));
        for j = 1:7         
            tmp(1:2:end,1:2:end,1:2:end) = Cup{k,j};
            Cup{k,j} =  circshift(imfilter(tmp,Psi(:,:,:,1),'corr','circular'),[-N/2+1,-N/2+1,-N/2+1]);
        end
    end
end
U = U+Cup{levels+1,1};