function y = enhancementTransform(x,a)

if nargin<2, a = 3; end
% y = 2./(1+exp(-a*x))-1;
% y = 2*sign(x)./(1 + exp(-a*abs(x))) - sign(x);
if a==0
    y = x;
else
    y = a*x./sqrt(1+(a*x).^2);
end
