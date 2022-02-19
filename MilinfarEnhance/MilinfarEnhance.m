function [U,out] = MilinfarEnhance(Z,M,opts)

[m,n,C,levels] = size(Z);
levels = levels-1;
if ~isfield(opts,'alpha')
    opts.alpha = 3;
end
if ~isfield(opts,'gamma')
    opts.gamma = .25;
end
if ~isfield(opts,'betas')
    opts.betas = [ones(1,levels-2)*2, 0.5, 0];
end


Zd = -diff(Z,1,4);
M = M.^opts.gamma;
U = Z(:,:,:,end);
out.Us = zeros(m,n,C,levels+1);
out.Us(:,:,:,1) = U;


cnt = 0;
for i = levels:-1:1
    cnt = cnt+1;
    tmp = enhancementTransform(Zd(:,:,:,i),opts.alpha);
    U = U + opts.betas(cnt)*tmp.*M;
    out.Us(:,:,:,cnt+1) = U;
end



