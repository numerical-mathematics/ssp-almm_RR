function [alpha,beta,tildebeta,C,tildeC] = ssp_lmm(k,p,method,downwind,y,solver)
% function [alpha,beta,tildebeta,C,tildeC] = ssp_lmm(k,p,method,downwind,y,solver)
%
% Author:  Yiannis Hadjimichael
% Created: November 2016
%%
% Finds optimal k-step SSP LMM methods with order of accuracy p that allow
% both upwind and/or downwind operators.
%
% The code works for classical or perturbed (i.e. downwind) LMMs and covers
% both explicit and implicit methods.
% 
% Input variables:
%   k                       -- number of steps
%   p                       -- order of accuracy
%   method                  -- explicit ('ex') or implicit ('im')
%   downwind                -- 'dw' for downwind, otherwise 'none'
%   y                       -- forward Euler step-size ratio hFE/tilde{hFE}
%                              hFE (tilde{hFE}) corresponds to the FE
%                              condition of upwind (downwind) operator 
%   solver                  -- optimization solver:
%                              'linprog', 'fmincon' or 'cvx'
%
% Output variables:
%   alpha, beta, tildebeta  -- method's coefficients
%   (C,tildeC)              -- SSP coefficients
%
% alpha is k x 1
% beta is k x 1 (explicit) or k+1 x 1 (implicit)
% tildebeta is k x 1 (explicit) or k+1 x 1 (implicit)
%
% Obtimization varaibles are stored in a single vector x:
% x = [gamma beta tildebeta],
% where gamma = alpha - (C*beta + tildeC*tildebeta) and
%
% Method is given by
%
% u_n = SUM_{j=0}^{k-1} alpha_j*u_{n-k+j} + 
%       h*SUM_{j=0}^{k}[beta_j*F(u_{n-k+j})-tildebeta_j*tildeF(u_{n-k+j})],
%
% where F(u) and tildeF(u) are upwind and downwind operators respectively.
%
% Notice 1:
% The 'active-set' algorithm used in this code with the 'linprog' solver
% will be replaced by 'interior-point' and 'dual-simplex' algorithm in
% future MATLAB releases. However, 'active-set' performs much better than
% the other two options, in particular for a large number of steps and
% order of accuracy. The relevant warning has been suppressed.
%
% Notice 2:
% For implicit methods with downwinding, 'linprog' solver may not provide
% the most accurate results, expecially when k>=12 and p>=10.
% This is because the matrix of the linear problem becomes ill-conditioned;
% in such cases it is advisable to use the nonlinear oprimization solver
% 'fmincon' again with the 'active-set' algorithm. This in general will
% result in more accurate results.
% By default the code uses the 'linprog' solver.
% More details about the specific cases that 'fmincon' gives better results
% can be found in 'check_SSP_coeff.m'. 

% =========================================================================

%% Sanity checks:
if strcmp(method, 'ex')
    implicit = 0;
elseif strcmp(method, 'im')
    implicit = 1;
else
    error('Choose ''ex'' for explicit or ''im'' for implicit.')
end
if strcmp(downwind, 'dw')
    downwind = 1;
elseif strcmp(downwind, 'none')
    downwind = 0;
else
    error('Choose ''dw'' for downind or ''none'' otherwise.')
end

% if ratio of FE step sizes is not provided, set it to unity and use
% 'linprog' as the defualt optimization solver
if nargin <= 5
    solver = 'linprog';
    if nargin == 4
        y = 1;
    end
end

% =========================================================================

%% Optimization routine:

ordertol = 1.e-11; % tolerance for satisfying order conditions
bisectol = 1.e-15; % tolerance for bisection termination

% optimization solver
if strcmp(solver,'linprog')
    algorithm = 'active-set';
    warning('off','optim:linprog:AlgOptsWillError')
elseif strcmp(solver,'fmincon')
    algorithm = 'active-set';
elseif ~strcmp(solver,'cvx')
    msg = ['Incorrect optimization solver. Choose one from: ' ...
        '''linprog'', ''fmincon'', ''cvx'''];
    error(msg)
end

% optimization options
if strcmp(solver,'cvx')
    cvx_precision high
    cvx_quiet true
else
    opts = optimoptions(solver,'Algorithm',algorithm,'MaxIter',1e6, ...
        'ConstraintTolerance',1.e-15,'Display','off');
end

% variables for vector of coefficients
%   without downwind:
%       explicit: 2*k
%       inplicit: 2*k+1
%   with downwind:
%       explicit: 3*k
%       implicit: 3*k+2
n = 2*k + implicit + downwind*(k + implicit);

Rmin = 0.0; % keep looking until a method with at least this value is found
Rmax = 2.01; % upper bound for radius of absolute monotonicity (r.a.m.)

% set the lower bounds on the unknowns
lb = zeros(1,n);
ub = zeros(1,n) + 1.e2;

% objective function
f = zeros(1,n);
ff = @(x) zeros(1,n)*x';
% ff = @(x) x*x';

% initial point for 'fmincon' solver
c = ones(1,n);

% preallocation
D = zeros(p+1,k); % for gamma
B = zeros(p+1,k+implicit); % for beta
if downwind
    tildeB = zeros(p+1,k+implicit); % for tildebeta
end
beq = ones(p+1,1); % RHS of order consitions
    
% bisection algorithm
rmin = Rmin;
rmax = Rmax;
r = rmax; % starting value for r

while (rmax - rmin > bisectol)
    
    % avoiding division by zero
    D(1,:) = 1.;
    B(1,1:end-implicit) = r;
    if downwind
        tildeB(1,1:end-implicit) = y*r;
    end
    
    % equality constraints (order conditions on gamma, beta, tildebeta)
    jk = linspace(0,k-1,k)/k;
    for i = 1:p
        jk_i = (jk).^i;
        jk_im1 = (jk).^(i-1);
        
        D(i+1,:) = jk_i;
        B(i+1,1:end-implicit) = r*jk_i + i/k*jk_im1;
        if downwind
            tildeB(i+1,1:end-implicit) = y*r*jk_i - i/k*jk_im1;
        end
        if implicit
            B(i+1,end) = i/k;
            if downwind
                tildeB(i+1,end) = -i/k;
            end
        end
    end
    
    if downwind
        Aeq = [D B tildeB];
    else
        Aeq = [D B];
    end
    
    % optimization program call
    if strcmp(solver,'linprog')
        [~,~,flag] = linprog(f,[],[],Aeq,beq,lb,ub,[],opts);
    elseif strcmp(solver,'fmincon')
        [~,~,flag] = fmincon(ff,c,[],[],Aeq,beq,lb,ub,[],opts);
    elseif strcmp(solver,'cvx')
        cvx_begin
            variable x(n)
            minimize(norm(x))
            subject to
                Aeq*x == beq
                0 <= x
        cvx_end
        if strcmp(cvx_status,'Solved')
            flag = 1;
        else
            flag = 0;
        end
    end
    
    % check feasibility
    if flag == 1
        rmin = r; % feasible solution found
    else
        rmax = r; % not feasible
    end
    
    % update r
    r = (rmax + rmin)/2;
    
end

% store current value of SSP coefficient
C = rmin;
tildeC = y*C;

%==========================================================================

%% Now need to find coefficients for optimal value of C:
% avoiding division by zero
D(1,:) = 1.;
B(1,1:end-implicit) = C;
if downwind
    tildeB(1,1:end-implicit) = tildeC;
end

% equality constraints (order conditions on gamma, beta, tildebeta)
for i = 1:p
    jk_i = (jk).^i;
    jk_im1 = (jk).^(i-1);
    
    D(i+1,:) = jk_i;
    B(i+1,1:end-implicit) = C*jk_i + i/k*jk_im1;
    if downwind
        tildeB(i+1,1:end-implicit) = y*C*jk_i - i/k*jk_im1;
    end
    if implicit
        B(i+1,end) = i/k;
        if downwind
            tildeB(i+1,end) = -i/k;
        end
    end
end

if downwind
    Aeq = [D B tildeB];
else
    Aeq = [D B];
end

% find coefficients for optimal method
if strcmp(solver,'linprog')
    [X,~,flag] = linprog(f,[],[],Aeq,beq,lb,ub,[],opts);
elseif strcmp(solver,'fmincon')
    [X,~,flag] = fmincon(ff,c,[],[],Aeq,beq,lb,ub,[],opts); X = X';
elseif strcmp(solver,'cvx')
    cvx_begin
        variable x(n)
        minimize(norm(x))
        subject to
            0.0 <= x
            Aeq*x == beq
    cvx_end
    cvx_status
    if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
        X = x;
        flag = 1;
    end
end

%==========================================================================

%% Output coefficients for optimal LMM with SSP coefficients (C,tildeC)
if flag ~= 1
    alpha = []; beta = []; tildebeta = [];
    msg = ['The optimization problem did not converge to a ' ...
        'feasible solution.'];
    warning(msg)
else
    gamma = X(1:k);
    beta = X(k+1:2*k+implicit);
    if downwind
        tildebeta = X(2*k+1+implicit:end);
        alpha = gamma + C*beta(1:end-implicit) + ...
            y*C*tildebeta(1:end-implicit);
        % obtain optimal method that satisfies beta'*tildebeta = 0
        % (see Lemma 2.7 in the paper
        % "Strong stability preserving additive linear multistep methods",
        % (Y. Hadjimichael and D. Ketcheson))
        delta = zeros(length(beta),1);
        tildedelta = zeros(length(tildebeta),1);
        for i = 1:length(beta)
            if beta(i)-tildebeta(i) > 0
                delta(i) = beta(i)-tildebeta(i);
            else
                tildedelta(i) = tildebeta(i)-beta(i);
            end
        end
        beta = delta;
        tildebeta = tildedelta;
    else
        tildebeta = [];
        alpha = gamma + C*beta(1:end-implicit);
    end
    
    % check if order conditions are satisfied
    order = check_order('lmm',alpha,beta,tildebeta,ordertol);
    if order ~= p
        msg = ['The optimal method does not satisfy the order ' ...
            'conditions. Consider increasing the tolerance (ordertol).'];
        warning(msg)
    end
end

end %function
