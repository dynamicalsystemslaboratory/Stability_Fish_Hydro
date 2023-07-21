%% Plot contour plots of natural frequencies of stable equilibria 
% in alpha-lambda plane

clear all; clc; close all
format short

FTsz = 25; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.1; rhostr = '0_1'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
tol = 1e-05; tol2=1e-3; 
tol4 = 0.1; % angle tol for IL
tol5 = 0.8; % angle tol for ST 
tol6 = 0.2; % angle tol for WF 
%%%%%%%%%%% When doing for rho=0 case, remove ialph=1 values!

% Save figure
save_fig_dir = './figures_natfreq/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
Fstr = '1Td_v6_fine2';
figname = sprintf('%s_%s_%s_%s_%s', 'Stable_w1_rho', rhostr, 'K', Kstr, Fstr);
figname2 = sprintf('%s_%s_%s_%s_%s', 'Stable_w2_rho', rhostr, 'K', Kstr, Fstr);

fdir = './files1Td_v6_fine2/';

fdir2 = [fdir 'om_fine/'];
fname = sprintf('%s%s_%s_%s_%s', fdir2, 'w1w2_rho', rhostr, 'K', Kstr);

% Read in parameters from one of the roots files
ilam = 1; ialph = 1; 
fdir1 = [fdir 'Roots3/'];
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
load(filename);

alpharr = [0.0:0.005:0.16 0.16+0.02:0.02:1.0]; 
lamarr = [0.0:0.01:1.5 1.5+0.02:0.02:3.0];
[Alph,Lam] = meshgrid(alpharr,lamarr);

filename = sprintf('%s%s',fname,'_B.mat');
load(filename);

F1=figure(1);
[c,h] = contourf(Alph,Lam,om1,50); 
colorbar; caxis([0 3.5]);
set(h, 'edgecolor','none');
xlabel('\alpha'); ylabel('\Lambda');
pbaspect([5 4 1]);
% title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
saveas(gcf,[save_fig_dir figname],'png');

% F2=figure(2);
% x= [0 1]; y = [0 3];
% % x = [0 1]; y = [0 1];
% imagesc(x,y,flipud(om1));
% set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
% colorbar; %caxis([0 4]);
% set(h, 'edgecolor','none');
% xlabel('\alpha'); ylabel('\Lambda');
% % title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
% saveas(gcf,[save_fig_dir figname '_imsc'],'png');

F3=figure(3);
[c,h] = contourf(Alph,Lam,om2,50); 
colorbar; %caxis([0 4]);
set(h, 'edgecolor','none');
xlabel('\alpha'); ylabel('\Lambda');
pbaspect([5 4 1]);
% title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
saveas(gcf,[save_fig_dir figname2],'png');

% F4=figure(4);
% x= [0 1]; y = [0 3];
% % x = [0 1]; y = [0 1];
% imagesc(x,y,flipud(om2));
% set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
% colorbar; %caxis([0 4]);
% set(h, 'edgecolor','none');
% xlabel('\alpha'); ylabel('\Lambda');
% % title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
% saveas(gcf,[save_fig_dir figname2 '_imsc'],'png');