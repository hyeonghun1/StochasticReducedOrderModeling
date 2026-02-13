function C = page_cov(x, transposePages)

% Compute covariance pagewise.
% x:  R^{n x L x s}
% xr: R^{r x L x s}


if transposePages
    x = pagetranspose(x);  % L x n x s
end

[~, n, s] = size(x);
C = zeros(n, n, s);

for ii=1:s
    C(:,:,ii) = cov(x(:,:,ii));
end

end




