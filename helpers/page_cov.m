function C = page_cov(x, transposePages)

% Compute covariance pagewise.


if transposePages
    x = pagetranspose(x);
end

[m, n, s] = size(x);
C = zeros(n, n, s);

for ii=1:s
    C(:,:,ii) = cov(x(:,:,ii));
end

end