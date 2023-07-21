%% This code plots the component of zero eigenvector of stable equilibria along lambda
clear all; clc; close all
format short

FTsz = 25; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

fdir = './files1Td_v6_fine2/';
filename2 = sprintf('%s%s',fdir, 'Ezero_fine.mat'); % To read

load(filename2)

save_fig_dir = './figures_natfreq/';

% Plot
figure;
[c,h] = contourf(Alph,Lam,Ezero); 
colorbar;
set(h, 'edgecolor','none');
xlabel('\alpha'); ylabel('\Lambda');
pbaspect([5 4 1]);
saveas(gcf,[save_fig_dir 'zeroEvec_lam'],'png');

% figure;
% x = [0 1]; y = [0 3];
% imagesc(x,y,flipud(Ezero));
% set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
% % colormap(cMap);
% colorbar; 
% % caxis([0 4]);
% xlabel('\alpha'); ylabel('\Lambda');
% pbaspect([5 4 1]);
% % title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
% saveas(gcf,[save_fig_dir 'zeroEvec_lam_imsc'],'png');


