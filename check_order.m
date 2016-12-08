function p = check_order(method,alpha,beta,tbeta,tol)
% function p = check_order(method,alpha,beta,tbeta,tol)
%
% Author:  Yiannis Hadjimichael
% Created: November 2016
%%
% Checks the order of accuracy for an explicit/implicit LMM or IMEX LMM
% given an error tolerance.
%
% Input variables:
%   method                  -- LMM ('lmm') or IMEX LMM ('imex')
%   alpha, beta, tildebeta  -- method's coefficients
%
% Output variables:
%   p                       -- order of accuracy

% =========================================================================

%%
% set default tolerance if not provided
if nargin == 4
    tol = 1.e-10;
end

k = length(alpha);
p = -1;
i = 1;
cond = abs(sum(alpha) - 1.); % consistency

if strcmp(method,'lmm')
    % for methods with downwind
    if ~isempty(tbeta)
        beta = beta - tbeta;
    end
    
    % check order conditions
    while cond < tol
        cond = 0.;
        for j=0:k-1
            cond = cond + alpha(j+1)*(j/k)^i + beta(j+1)*i/k*(j/k)^(i-1);
        end
        cond = cond - 1.;
        if length(beta) == k+1
            cond = cond + beta(end)*i/k;
        end
        cond = abs(cond);
        
        p = p + 1;
        i = i + 1;
    end
elseif strcmp(method,'imex')
    % check order conditions
    cond2 = 0.;
    while cond < tol && cond2 < tol
        cond = 0.;
        cond2 = 0.;
        for j=0:k-1
            cond = cond + alpha(j+1)*(j/k)^i + beta(j+1)*i/k*(j/k)^(i-1);
            cond2 = cond2 + alpha(j+1)*(j/k)^i + tbeta(j+1)*i/k*(j/k)^(i-1);
        end
        cond = cond - 1.;
        cond2 = cond2 + tbeta(end)*i/k - 1.;
        cond = abs(cond);
        cond2 = abs(cond2);
        
        p = p + 1;
        i = i + 1;
    end
end

if p == -1
    warning('The method is not consistent for the given tolerance.')
end

end %function

