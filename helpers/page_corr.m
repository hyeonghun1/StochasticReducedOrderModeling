function R = page_corr(x, transposePages)
% Compute correlation pagewise from 3D array
% x:  n x L x s (or L x n x s if transposePages = true)
% R:  n x n x s, correlation matrices

if transposePages
    x = pagetranspose(x);  % L x n x s
end

[~, n, s] = size(x);
R = zeros(n, n, s);

for ii = 1:s
    C = cov(x(:,:,ii));         % covariance matrix
    D = diag(1 ./ sqrt(diag(C))); % normalization matrix
    R(:,:,ii) = D * C * D;      % correlation matrix
end

end
