function [E,A,B,N,M,f,n,m,h,s, ind_solution] = heat2dExample()
    n = 34;
    dx = 1/(n+1);
    h = 1e-2;
    mu1 = 1e-2;
    mu2 = 0;      % 1e-1;
    s = 100;
  
    % 1d Laplace
    A1d = mu1*sparse(fd_laplace(n,dx));
    %A1d = A1d + mu2*fd_advection(n,dx);
  
    % 2d Laplace
    Afull = kron(speye(n),A1d)+kron(A1d,speye(n));

    % define coords of the rectangle that is removed
    cwidth = 10;
    clength = 14;
    cstart = n/2;
    cind = ([cstart:cstart+cwidth-1]) + ([cstart:cstart+clength-1]'*n);

    % solution is only computed on these indices. the complement is 0.
    ind_solution = setdiff(1:n^2, cind);
   
    % dynamics in rectangle removed
    A = Afull(ind_solution,ind_solution);
    E = speye(size(A));

    % define input source rectangle
    istart = 1;
    k = 12;
    Bfull = zeros(n^2, 1);
    for ii = cstart - k : cstart + k
        ind_source = (istart:(istart+k-1)) + n*(ii-1);
        Bfull(ind_source) = 1;
    end
    
    B = sparse(Bfull(ind_solution,:));
    [n, m] = size(B);

    % no bilinearity
    N = sparse(zeros(size(A)));

    % the 1d noise acts on the same rectangle as the input
    % M = 0.01*B;
    M = B/norm(B);
    
	
    %M = zeros(size(M));
    K = 1;
  
    % drift implicit Euler-Maruyama
    % f = @(x0,u,L) (E - h*A - h*N*u)\(x0 + h*B*u ...
    %                                     + sqrt(h)*M*randn(size(M,2), L));
    f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h);
end

function x_next = EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)
    
    % Deterministic term
    Bu = B * u;
    det_term = x0 + h*Bu;

    % Stochastic term
    xi = randn(size(M, 2), L);
    rand_term = sqrt(h) * M * xi;

    % RHS
    rhs = det_term + rand_term;

    % Solve linear system
    x_next = (E - h*A - h*N*u) \ rhs;
end
