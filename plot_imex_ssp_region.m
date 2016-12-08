function [C,hatC] = plot_imex_ssp_region(k,p)
% function [C,hatC] = plot_imex_ssp_region(k,p)
%
% Author:  Yiannis Hadjimichael
% Created: November 2016
%%
% Computes the feasibility region of IMEX SSP LMMs for a given number of
% steps k and order of accuracy p.
%
% Input variables:
%   k                       -- number of steps
%   p                       -- order of accuracy
%
% Output variables:
%   (C,hatC)                -- SSP coefficients
%   plot of SSP region

% =========================================================================

%%
% ratio of forward Euler step sizes of explicit and implicit parts
y = tan(0:pi/100:pi/2);
y(end) = tan(pi/2-1.e-3);

% compute optimal SSP coefficients for each value of y
C = zeros(length(y),1);
hatC = zeros(length(y),1);
for i = 1:length(y)
    [~,~,~,C(i),hatC(i)] = ssp_imex(k,p,y(i));
end

% plotting
fig = figure();
plot(C,hatC,'-','LineWidth',2);
hc = get(fig,'children'); set(hc,'fontsize',18);
xlabel('$$r$$','FontSize',24,'Interpreter','latex');
y = ylabel('$$\hat{r}$$','FontSize',24,'Interpreter','latex','Rotation',0);
pos = get(y,'position');
set(y,'Units','Normalized','Position',[-(0.16+pos(1)) 0.45 -1]);
get(y,'position')
grid on;

% saving figure
file = sprintf('/figures/IMEX(%d,%d).pdf',k,p);
saveas(fig,[pwd file]);

end %function
