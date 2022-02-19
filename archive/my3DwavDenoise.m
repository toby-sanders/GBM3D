function X = my3DwavDenoise(X,tau,opts)

% basic 3D wavelet denoising
if nargin<3
    opts.wname = 'sym4';
    opts.levels = 4;
else
    if ~isfield(opts,'wname'), opts.wname = 'sym4'; end
    if ~isfield(opts,'levels'), opts.levels = 4; end
end

X = wavedec3(X,opts.levels,opts.wname); % decompose
for i = 1:numel(X.dec) % threshold coefficients
    X.dec{i} = X.dec{i}.*(abs(X.dec{i})>tau); 
end
X = waverec3(X); % reconstruct