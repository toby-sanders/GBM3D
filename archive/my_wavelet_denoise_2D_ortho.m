function U = my_wavelet_denoise_2D_ortho(wname,levels,U,tau)

% Written by Toby Sanders
% orthonormal wavelet decomposition and reconstruction

% these are designed for 1-D problems

% inputs:
% wname - name of the wavelet used
% levels - number of levels in the decomposition
% U - signal to be denoised
% tau - thresholding value


[p,q,r] = size(U);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);
Lo_d2 = Lo_d'*Lo_d;
Hi_d2 = Hi_d'*Hi_d;
HL1 = Hi_d'*Lo_d;
HL2 = Lo_d'*Hi_d;

figure(134);
subplot(2,2,1);imagesc(Lo_d2);
subplot(2,2,2);imagesc(Hi_d2);
subplot(2,2,3);imagesc(HL1);
subplot(2,2,4);imagesc(HL2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet decomposition step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = cell(levels+1,3); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    for j = 1:3
        if j==1
            C{i,j} = circshift(imfilter(U,Hi_d2,'conv'),[N/2,N/2]);
        elseif j==2
            C{i,j} = circshift(imfilter(U,HL1,'conv'),[N/2,N/2]);
        else
            C{i,j} = circshift(imfilter(U,HL2,'conv'),[N/2,N/2]);
        end
        C{i,j} = C{i,j}(1:2:end,1:2:end);
    end
    
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,Lo_d2,'conv'),[N/2,N/2]);
    U = U(1:2:end,1:2:end);
end
C{levels+1,1} = U; % last step returns approximation coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% throw away small detail coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:levels
    for j = 1:3
        C{i,j} = max(abs(C{i,j})-tau,0).*sign(C{i,j});
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet reconstruction algorithm, which performs the adjoint operations
% to the above loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(p,q); % store restored signal into U


for i = 1:levels
    for k = 1:3
       tmp = zeros(size(C{i,k},1)*2,size(C{i,k},2)*2);
       tmp(1:2:end,1:2:end) = C{i,k};
       if k==1
           tmp = circshift(imfilter(tmp,Hi_d2,'corr'),[-N/2+1,-N/2+1]);
       elseif k==2
           tmp = circshift(imfilter(tmp,HL1,'corr'),[-N/2+1,-N/2+1]);
       else
           tmp = circshift(imfilter(tmp,HL2,'corr'),[-N/2+1,-N/2+1]);
       end
       for j = 1:i-1
           tmp2 = zeros(size(tmp,1)*2,size(tmp,2)*2);
           tmp2(1:2:end,1:2:end) = tmp;
           tmp =  circshift(imfilter(tmp2,Lo_d2,'corr'),[-N/2+1,-N/2+1]);
       end
    
       U = U + tmp(1:p,1:q);
    end
end

tmp = C{levels+1,1};
for i = 1:levels
    tmp2 = zeros(size(tmp,1)*2,size(tmp,2)*2);
    tmp2(1:2:end,1:2:end) = tmp;
    tmp = circshift(imfilter(tmp2,Lo_d2,'corr'),[-N/2+1,-N/2+1]);
end
U = U + tmp(1:p,1:q);


