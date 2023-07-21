%% This code plots heatmap of In-line (IL), Staggered (ST), Wall (W) stable
%% or all unstable (or other soln)
clear all; clc; close all
format short

FTsz = 20; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.0; rhostr = '0'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
tol = 1e-05; tol2=1e-3; 
tol4 = 0.1; % angle tol for IL
tol5 = 0.8; % angle tol for ST 
tol6 = 0.2; % angle tol for WF 

%%%%%%%%%%% When doing for rho=0 case, remove ialph=1 values!

%% Set up alph, lam arrays - non-uniform
alpharr = [0.0:0.005:0.16 0.16+0.02:0.02:1.0]; 
lamarr = [0.0:0.01:1.5 1.5+0.02:0.02:3.0];
[Alph,Lam] = meshgrid(alpharr,lamarr);
flagmesh = zeros(size(Alph));
Nalph1 = size(0.0:0.005:0.16,2); Nalph2 = size(0.16+0.02:0.02:1.0,2);
Nlam1 = size(0.0:0.01:1.5,2); Nlam2 = size(1.5+0.02:0.02:3.0,2);
Nlam = Nlam1+Nlam2; Nalph = Nalph1+Nalph2;
flagmesh(1:Nlam1,1:Nalph1) = 1; % Folder Roots3
flagmesh(Nlam1+1:Nlam,1:Nalph1) = 2; % Folder Roots3b
flagmesh(1:Nlam1,Nalph1+1:Nalph) = 3; % Folder Roots3c

% Folder Roots
lamarr_coarse = 0.0:0.02:3.0; Nlam_coarse = size(lamarr_coarse,2);
alpharr_coarse = 0.0:0.02:1.0; Nalph_coarse = size(alpharr_coarse,2);
Nlam1_coarse = size(0:0.02:1.5,2); Nalph1_coarse = size(0:0.02:0.16,2);

fdir = './files5Td_v6_fine2/';
fdirB = './../MATLAB_CODES/files5Td_v6_fine2/';

% Read in parameters from one of the roots files
ilam = 1; ialph = 1; 
fdir1 = [fdir 'Roots3/'];
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
load(filename);

% Save figure
save_fig_dir = './figures/';
Fstr = '5Td_v6_fine2';
figname = sprintf('%s_%s_%s_%s_%s', 'Stab_rho', rhostr, 'K', Kstr, Fstr);

%% ilam, ialph: global indices; Local indices to read from file are different
for ilam = 1:Nlam
    if (r==0)
        ialph=1; N_IL_MS(ilam,ialph) = 1; ialphst = 2;
    else
        ialphst = 1;
    end

    for ialph = ialphst:Nalph
        l = Lam(ilam,ialph); lam = l;
        a = Alph(ilam,ialph);
%         disp(a)

        if (flagmesh(ilam,ialph)==1)
            fdir1 = [fdir 'Roots3/'];
            fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
            filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
        elseif (flagmesh(ilam,ialph)==2)
            fdir1 = [fdir 'Roots3b/'];
            fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
            filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam-Nlam1, 'ialph', ialph, '.mat');
        elseif (flagmesh(ilam,ialph)==3)
            fdir1 = [fdir 'Roots3c/'];
            fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
            filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph-Nalph1, '.mat');
        else
            fdir1 = [fdirB 'Roots/'];
            ilam_coarse = find(abs(lamarr_coarse-l)<1e-5); ialph_coarse = find(abs(alpharr_coarse-a)<1e-5);
            fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
            filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam_coarse, 'ialph', ialph_coarse, '.mat');
        end
        load(filename);

        if (isempty(Roots3))
            continue;
        else
            eigflag = Roots3(:,6);

            %%%% eigflag = -1:U, 0:MS, 1:S 
            % Stable (S)
            idxS = eigflag > 0; N_S(ilam,ialph) = sum(idxS);
            
            % Marg. Stable (MS)
            idxMS = eigflag == 0; N_MS(ilam,ialph) = sum(idxMS);
            
            % Unstable (U)
            idxU = eigflag < 0; N_U(ilam,ialph) = sum(idxU);

            % Ups or Downstream
            idxUpD = ( (abs(Roots3(:,3)) < tol4) & (abs(Roots3(:,4)) < tol4) ) | ...
                     ( (abs(abs(Roots3(:,3))-pi) < tol4) & (abs(abs(Roots3(:,4))-pi) < tol4) );
            idxUpD2 = ( (abs(Roots3(:,3)) < tol5) & (abs(Roots3(:,4)) < tol5) ) | ...
                     ( (abs(abs(Roots3(:,3))-pi) < tol5) & (abs(abs(Roots3(:,4))-pi) < tol5) );
            
            % Wall facing
            idxWf =  ( (abs(abs(Roots3(:,3))-pi/2) < tol6) & (abs(abs(Roots3(:,4))-pi/2) < tol6) );

            % IL
            idxIL = abs(Roots3(:,1)-Roots3(:,2))<tol2;
            idx_IL_S = (idxS==1 & idxIL==1 & idxUpD==1); N_IL_S(ilam,ialph) = sum(idx_IL_S);
            idx_IL_MS = (idxMS==1 & idxIL==1 & idxUpD==1); N_IL_MS(ilam,ialph) = sum(idx_IL_MS);
            
            % ST (staggered symm about centerline)
            idxST = ( abs(Roots3(:,1)+Roots3(:,2))<tol2 ) & abs(Roots3(:,1))>tol2 & abs(Roots3(:,2))>tol2;
            idx_ST_S = (idxS==1 & idxST==1 & idxUpD2==1); N_ST_S(ilam,ialph) = sum(idx_ST_S);
            idx_ST_MS = (idxMS==1 & idxST==1 & idxUpD2==1); N_ST_MS(ilam,ialph) = sum(idx_ST_MS);
            
            % Other stable config
            idx_W_S = (idxS==1 & idxWf==1); N_W_S(ilam,ialph) = sum(idx_W_S);
            idx_W_MS = (idxMS==1 & idxWf==1); N_W_MS(ilam,ialph) = sum(idx_W_MS);

            % count total S eq pts
            countcat = (sum(idx_IL_S)+sum(idx_IL_MS)+sum(idx_ST_S)+sum(idx_ST_MS)+sum(idx_W_S)+sum(idx_W_MS));
            count = sum(idxS) + sum(idxMS);

            if (countcat<count)
                disp('Error! Some S eq pts are not being counted ...');
%                 Roots3
%                 pause;
            end
        end

    end
    disp(['Done lambda = ' num2str(l)]);
end


%% Plot
flagid = zeros(size(Alph));
for ilam = 1:size(Alph,1)
    for ialph = 1:size(Alph,2)
        l = Lam(ilam,ialph); lam = l;
        a = Alph(ilam,ialph);
        if ( (N_IL_MS(ilam,ialph)>0 || N_IL_S(ilam,ialph)>0) && ...
                (N_W_MS(ilam,ialph)>0 || N_W_S(ilam,ialph)>0) )
            flagid(ilam,ialph) = 4; % Both IL and W
        elseif (N_IL_MS(ilam,ialph)>0 || N_IL_S(ilam,ialph)>0)
            flagid(ilam,ialph) = 1; % IL only
        elseif (N_ST_MS(ilam,ialph)>0 || N_ST_S(ilam,ialph)>0)
            flagid(ilam,ialph) = 2; % ST only
        elseif (N_W_MS(ilam,ialph)>0 || N_W_S(ilam,ialph)>0)
            flagid(ilam,ialph) = 3; % W only
        else
            flagid(ilam,ialph) = 0; % None stable
        end
    end
end

%colors:
cMap = [0, 0, 0; ... % White for 0 (none)
        0, 0, 1; ... % Blue  for 1 (IL)
        0, 1, 1; ... % ... for 2 (ST)
        0, 1, 0; ... % Green for 3 (W)
        1, 0, 0];    % Red   for 4 (Both IL-W)

F1=figure(1);
[c,h] = contourf(Alph,Lam,flagid); 
colormap(cMap);
cb = colorbar; caxis([0 4]);
set(h, 'edgecolor','none');
xlabel('\alpha'); ylabel('\Lambda');
cb.Ticks = linspace(0,4,5); 
cb.TickLabels{1} = 'None';
cb.TickLabels{2} = 'IL';
cb.TickLabels{3} = 'ST';
cb.TickLabels{4} = 'W';
cb.TickLabels{5} = 'IL,W';
title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
saveas(gcf,[save_fig_dir figname],'png');

% figure;
% x = [0 1]; y = [0 3];
% imagesc(x,y,flipud(flagid));
% set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
% colormap(cMap);
% cb = colorbar; caxis([0 4]);
% xlabel('\alpha'); ylabel('\Lambda');
% cb.Ticks = linspace(0,4,5); 
% cb.TickLabels{1} = 'None';
% cb.TickLabels{2} = 'IL';
% cb.TickLabels{3} = 'ST';
% cb.TickLabels{4} = 'W';
% cb.TickLabels{5} = 'IL,W';
% title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
% saveas(gcf,[save_fig_dir figname '_imsc'],'png');

