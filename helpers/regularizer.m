function Reg = regularizer(r, lambda, isbilinear)

%REGULARIZER Construct regularizer diagonal matrix
%
% Inputs:
%   r        - reduced dimension (integer)
%   lambda1 - linear regularizer (scalar)
%
% Output:
%   Reg     - diagonal regularizer matrix

% d = 1 + r + r*(r+1)/2;     % total dimension

m = 1;

if isbilinear
    d = r + m + r*m;
else
    d = r + m;
end

Reg_vec = zeros(d, 1);
Reg_vec(1:r) = lambda(1);   % Regularizer for A
Reg_vec(r+1) = lambda(2);   % Regularizer for B
if isbilinear
    Reg_vec(r+2:end) = lambda(3);  % Regularizer for N
end

Reg = diag(Reg_vec);

end