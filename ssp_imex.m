function [alpha,beta,hatbeta,C,hatC] = ssp_imex(k,p,y)
% function [alpha,beta,hatbeta,C,hatC] = ssp_imex(k,p,y)
%
% Author:  Yiannis Hadjimichael
% Created: November 2016
%%
% Finds optimal k-step IMEX methods with order of accuracy p for which
% both explicit and implicit parts are SSP.
% 
% Input variables:
%   k                       -- number of steps
%   p                       -- order of accuracy
%   y                       -- forward Euler step-size ratio hFE/hat{hFE}
%                              hFE (hat{hFE}) corresponds to the FE
%                              condition of nonstiff (stiff) operator 
%
% Output variables:
%   alpha, beta, tildebeta  -- method's coefficients
%   (C,hatC)                -- SSP coefficients
%
% alpha is k x 1
% beta is k x 1 (explicit part)
% hatbeta is k+1 x 1 (implicit part)
%
% Obtimization variables are stored in a single vector x:
% x = [gamma beta hatbeta].
% Here gamma = alpha - (r*beta + hatr*hatbeta) and
% r <= C, hatr <= hatC, where (C,hatC) are the SSP coefficients.
%
% Method is given by
%
% u_n = SUM_{j=0}^{k-1} alpha_j*u_{n-k+j} + 
%       h*SUM_{j=0}^{k-1} beta_j*F(u_{n-k+j}) + 
%       h*SUM_{j=0}^{k} hatbeta_j*hatF(u_{n-k+j})],
%
% where F(u) and hatF(u) are the operators solved explicitly and implicitly
% respectively.
%
% Notice 1:
% The 'active-set' algorithm used in this code with the 'linprog' solver
% will be replaced by 'interior-point' and 'dual-simplex' algorithm in
% future MATLAB releases. However, 'active-set' performs much better than
% the other two options, in particular for a large number of steps and
% order of accuracy. The relevant warning has been suppressed.

% =========================================================================

%% Editable options:
ordertol = 1.e-12; % tolerance for satisfying order conditions
bisectol = 1.e-15; % tolerance for bisection termination

% if ratio of FE step sizes is not provided, set it to unity by default
if nargin == 2
    y = 1;
end

%==========================================================================

%% Optimization routine

% set optimization parameters:
opts = optimoptions(@linprog,'Algorithm','active-set','MaxIter',1e6, ...
    'ConstraintTolerance',1.e-15,'Display','none');
warning('off','optim:linprog:AlgOptsWillError')% suppress algorithm warning
           
% variables for vector of coefficients
n = 3*k + 1;

Rmin = 0.0; % keep looking until a method with at least this value is found
Rmax = n+1.e-2; % upper bound for radius of absolute monotonicity (r.a.m.)

% set the lower bounds on the unknowns
lb = zeros(1,n);

% dumping objective function
f = zeros(1,n);

% preallocation
D = zeros(2*p+1,k); % for gamma's
B = zeros(2*p+1,k); % for beta's
hatB = zeros(2*p+1,k+1); % for betahat's
beq = zeros(2*p+1,1); % RHS of order consitions
    
% bisection algorithm
rmin = Rmin;
rmax = Rmax;
r = rmax; % starting value for r
while (rmax - rmin > bisectol)
    % avoiding division by zero
    D(1,:) = 1.;
    D(p+1:end,:) = 0.;
    B(1,1:end) = r;
    hatB(1,1:end-1) = y*r;
    hatB(p+2:end,end) = -1;
    beq(1:p+1) = 1.;
    
    % equality constraints (order conditions on gamma, beta, hatbeta)
    jk = linspace(0,k-1,k)/k;
    for i=1:p
        jk_i = (jk).^i;
        jk_im1 = (jk).^(i-1);
        
        D(i+1,:) = jk_i;
        B(i+1,:) = r*jk_i + i/k*jk_im1;
        B(p+1+i,:) = jk_im1;
        hatB(i+1,1:end-1) = y*r*jk_i;
        hatB(p+1+i,1:end-1) = -jk_im1;
    end
    
    Aeq = [D B hatB];
    
    % the optimization call
    [~,~,flag] = linprog(f,[],[],Aeq,beq,lb,[],[],opts);
    
    % check feasibility
    if flag == -2
        rmax = r; % not feasible
    else
        rmin = r; % feasible solution found
    end
    
    % update r
    r = (rmax + rmin)/2;
    
end

% store current value of SSP coefficient
C = rmin;
hatC = y*C;

%==========================================================================

%% Now need to find coefficients for optimal value of C:
% avoiding division by zero
D(1,:) = 1.;
D(p+1:end,:) = 0.;
B(1,1:end) = C;
hatB(1,1:end-1) = hatC;
hatB(p+2:end,end) = -1;
beq(1:p+1) = 1.;

% equality constraints (order conditions on gamma, beta, hatbeta)
jk = linspace(0,k-1,k)/k;
for i=1:p
    jk_i = (jk).^i;
    jk_im1 = (jk).^(i-1);
    
    D(i+1,:) = jk_i;
    B(i+1,:) = C*jk_i + i/k*jk_im1;
    B(p+1+i,:) = jk_im1;
    hatB(i+1,1:end-1) = hatC*jk_i;
    hatB(p+1+i,1:end-1) = -jk_im1;
end

Aeq = [D B hatB];

% find coefficients for optimal method
[X,~,flag] = linprog(f,[],[],Aeq,beq,lb,[],[],opts);

%==========================================================================

%% Output coefficients for optimal IMEX LMM with SSP coefficients (C,hatC)
if flag ~= 1
    alpha = []; beta = []; hatbeta = [];
    msg = ['The optimization problem did not converge to a ' ...
        'feasible solution.'];
    warning(msg)
else
    gamma = X(1:k);
    beta = X(k+1:2*k);
    hatbeta = X(2*k+1:end);
    alpha = gamma + C*beta(1:end) + hatC*hatbeta(1:end-1);
    
    % check if order conditions are satisfied
    order = check_order('imex',alpha,beta,hatbeta,ordertol);
    if order ~= p
        msg = ['The optimal method does not satisfy the order ' ...
            'conditions. Consider increasing the tolerance (ordertol).'];
        warning(msg)
    end
end

end %function
