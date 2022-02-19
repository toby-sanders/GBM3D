function tau = getColoredTau(levels,sigmaPSD,tau0)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

% these thresholds are for colored denosing
% sigmaPSD = sigma^2 * PSD, where the power spectrum is the squared 
% magnitude of the Fourier transform of the filter

tau = zeros(levels,7);

% constants needed for wavelet filters
Ca = 0.016572815184060;
Cb = 0.121533978016438;
Cc = sqrt(2)/2; 

% bi-orthogonal decomposition wavelet filters
Lo_dxy = [Ca,-Ca,-Cb,Cb,Cc,Cc,Cb,-Cb,-Ca,Ca];
Hi_dxy = zeros(1,10); 
Hi_dxy(5) = -Cc; Hi_dxy(6) = Cc;

% build all 8 orthogonal 2D filters from the 1D filters
Psi = (zeros(10,10,8,'single'));
Psi(:,:,1) =  Lo_dxy'*Lo_dxy;
Psi(:,:,2) =  Hi_dxy'*Lo_dxy;
Psi(:,:,3) =  Lo_dxy'*Hi_dxy;
Psi(:,:,4) =  Lo_dxy'*Lo_dxy;
Psi(:,:,5) =  Hi_dxy'*Hi_dxy;
Psi(:,:,6) =  Hi_dxy'*Lo_dxy;
Psi(:,:,7) =  Lo_dxy'*Hi_dxy;
Psi(:,:,8) =  Hi_dxy'*Hi_dxy;

% Psi(:,:,1) =  Lo_dxy'*Lo_dxy;
% Psi(:,:,2) =  Hi_dxy'*Lo_dxy;
% Psi(:,:,3) =  Lo_dxy'*Hi_dxy;
% Psi(:,:,4) =  Lo_dxy'*Hi_dxy;
% Psi(:,:,5) =  Hi_dxy'*Hi_dxy;
% Psi(:,:,6) =  Hi_dxy'*Hi_dxy;
% Psi(:,:,7) =  Hi_dxy'*Hi_dxy;
% Psi(:,:,8) =  Hi_dxy'*Hi_dxy;


[d1,d2] = size(sigmaPSD);
for i = 1:levels
    for j = 1:7
        % compute normalized inner product of PSD with the wavelet filters
        PsiHat = imresize(Psi(:,:,j+1),2^(i-1),'nearest');
        PsiHat = abs(fft2(PsiHat,d1,d2)).^2;
        tau(i,j) = sqrt(sum(col(sigmaPSD.*PsiHat))./sum(col(PsiHat)));
    end
end
% empirically found this constant scaling factor is suitable
if nargin<3
    tau = tau*6;
else
    tau = tau*tau0;
end