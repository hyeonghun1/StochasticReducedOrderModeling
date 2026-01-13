function Reg = regularizer(r, lambda)

%REGULARIZER Construct regularizer diagonal matrix
%
% Inputs:
%   r        - reduced dimension (integer)
%   lambda1 - linear regularizer (scalar)
%
% Output:
%   Reg     - diagonal regularizer matrix

% d = 1 + r + r*(r+1)/2;     % total dimension

d = r;

Reg_vec = zeros(d, 1);
Reg_vec(1:r) = lambda; % Regularizer for A
% Reg_vec(r+2:end) = lambda2; % Regularizer for H

Reg = diag(Reg_vec);

end