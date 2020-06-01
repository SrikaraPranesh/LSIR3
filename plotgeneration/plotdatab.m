
function plotdatab(n1,n2,svnm,kappa,u,showleg,xmax)

% PLOTDATAAB creates a plot for the two experiments whose output data is 
% stored in MATLAB data objects n1 and n2. 
%
% Input
%       n1:         .mat file name for output of LSIR run
%       n2:         .mat file name for output of G-LSIR run
%       svnm:       filename for saving the generated plot
%       kappa:      value of 2-norm condition number of problem (for title of plot)
%       u:          value of working precision u for this problem
%       showleg:    whether to display legend (0 or 1)
%       xmax:       maximum x-axis limit


% Load data; extract xerr and rerr vectors
S = load(strcat('data/',n1,'.mat'));
G = load(strcat('data/',n2,'.mat'));
SNS = fieldnames(S);
SNG = fieldnames(G);
xerr1 = S.(SNS{1}).xerr;
rerr1 = S.(SNS{1}).rerr;
xerr2 = G.(SNG{1}).xerr;
rerr2 = G.(SNG{1}).rerr;

% Generate plot
f = figure();
semilogy(0:numel(xerr1)-1, xerr1, 'bo-','MarkerSize',10,'LineWidth',1);
hold on
semilogy(0:numel(rerr1)-1, rerr1, 'bx-','MarkerSize',10,'LineWidth',1);
hold on
semilogy(0:numel(xerr2)-1, xerr2, 'ro--','MarkerSize',10,'LineWidth',1);
hold on
semilogy(0:numel(rerr2)-1, rerr2, 'rx--','MarkerSize',10,'LineWidth',1);
hold on 
semilogy(0:xmax, u.*ones(xmax+1,1), 'k--','MarkerSize',10,'LineWidth',1);

xlabel('refinement step','FontSize',16,'Interpreter','latex')
if (showleg)
legend({'LSIR $x$', 'LSIR $r$', 'G-LSIR $x$', 'G-LSIR $r$'},'Interpreter','latex')
end
title(strcat('$\kappa_{\infty}(A) = $ ', num2str(kappa, '%1.e\n')),'FontSize',16,'Interpreter','latex')

axis([0, xmax, 1e-2*u, 10*double(max(max(xerr1),max(xerr2)))])
set(gca,'FontSize',16)

% Ensure only integers labeled on x axis
    atm = get(gca,'xticklabels');
    m = str2double(atm);
    xlab = [];
    num = 1;
    for i = 1:numel(m)
        if ceil(m(i)) == m(i)
            xlab(num) = m(i);
            num = num + 1;
        end
    end
    set(gca,'xticklabels',xlab);
    set(gca,'xtick',xlab);

% Save plot as fig and as pdf
saveas(f, strcat('fig/',svnm), 'fig');
saveas(f, strcat('fig/',svnm), 'pdf');

end