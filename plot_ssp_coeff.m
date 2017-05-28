function C = plot_ssp_coeff(k,p,method,solver)
% function C = plot_ssp_coeff(k,p,method,solver)
%
% Author:  Yiannis Hadjimichael
% Created: November 2016
%%
% Plot the SSP coefficients curves w.r.t the forward Euler step-size ratio
% hFE/tilde{hFE} for a given number of steps k and order of accuracy p.
%
% Input variables:
%   k                       -- number of steps
%   p                       -- order of accuracy
%   method                  -- explicit ('ex') or implicit ('im')

% Output variables:
%   C                       -- SSP coefficient C(y)
%   plot of SSP coefficients w.r.t y = hFE/tilde{hFE}

% =========================================================================

%%
if nargin == 3
    solver = 'linprog';
end

% ratio of forward Euler step sizes of upwind and downwind operators
y = linspace(0,1.5,30);

% compute optimal SSP coefficients for each value of y
C = zeros(length(y),1);
tildeC = zeros(length(y),1);
for i = 1:length(y)
    [~,~,~,C(i),tildeC(i)] = ssp_lmm(k,p,method,'dw',y(i),solver);
end
[~,~,~,R,~] = ssp_lmm(k,p,method,'dw',1,solver);

% plotting
fig = figure();
plot(y,C,y,tildeC,'-.','LineWidth',2);
hold on;
plot(y, R*ones(length(y)),'--k','LineWidth',2);
ssp = sprintf('$\\mathcal{C}_{%d,%d}(\\xi)$',k,p);
tildessp = sprintf('$\\tilde{\\mathcal{C}}_{%d,%d}(\\xi)$',k,p);
hc = get(fig,'children'); set(hc,'fontsize',18);
hleg = legend({ssp,tildessp},'Interpreter','Latex', ...
    'FontSize',20,'Location','SouthEast');
leg_pos = get(hleg,'position');
leg_pos = [leg_pos(1)-0.1*leg_pos(3) leg_pos(2)+0.1*leg_pos(4) ...
    1.1*leg_pos(3) 1.5*leg_pos(4)];
set(hleg,'position',leg_pos);
xlabel('$\xi$','FontSize',22,'Interpreter','Latex');
grid on;

% saving figure
file = sprintf('/figures/%s_DLMM(%d,%d).pdf',method,k,p);
saveas(fig,[pwd file]);

end %function
