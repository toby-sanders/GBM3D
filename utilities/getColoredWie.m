function sigmaWie = getColoredWie(sigmaPSD,d)

% Written by Toby Sanders @Lickenbrock Tech
% Date: 5-7-21

[~,ghat] = makeGausPSF([size(sigmaPSD,1),size(sigmaPSD,2)],3);
sigmaPSD = real(ifft2(ghat.*fft2(sigmaPSD)));
sigmaPSD = imresize(sigmaPSD,[d,d],'method','bilinear');
sigmaWie = zeros(d,d);
x = zeros(d);
for i = 1:d
    for j = 1:d
        x(:) = 0;
        x(j,i) = 1;
        y = idct2(x);
        y = abs(fft2(y)).^2;
        sigmaWie(j,i) = sqrt(sum(sigmaPSD(:).*y(:)))/d;
    end
end