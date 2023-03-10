%% Plot heatmap if In-line (IL), Staggered (ST), IL-ST stable, 
% or all unstable (or other soln)
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
r=0.1; rhostr = '0_1'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
tol = 1e-05; tol2=1e-3; 
tol4 = 0.1; % angle tol for IL
tol5 = 0.8; % angle tol for ST 
tol6 = 0.2; % angle tol for WF 

%%%%%%%%%%% When doing for rho=0 case, remove ialph=1 values!

% Read in parameters from one of the roots files
ilam = 1; ialph = 1; 
fdir1 = './files1Td_v6_fine2/Roots/';
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
load(filename);

lamarr = 0.0:0.02:3.0; Nlam = size(lamarr,2);
alpharr = 0.0:0.02:1.0; Nalph = size(alpharr,2);
[Alph,Lam] = meshgrid(alpharr,lamarr);

% Save figure
save_fig_dir = './figures6/';
Fstr = '1Td_v6_fine2';
figname = sprintf('%s_%s_%s_%s_%s', 'Stab_rho', rhostr, 'K', Kstr, Fstr);
% figname = sprintf('%s_%s_%s_%s', 'Num_CorrNoWall_rho', rhostr, 'K', Kstr);
% figname = sprintf('%s_%s_%s_%s%s', 'Num_Corr_rho', rhostr, 'K', Kstr, '_IC2');
% figname = sprintf('%s_%s_%s_%s%s', 'Num_CorrNoWall_rho', rhostr, 'K', Kstr, '_IC2');

%%

for ilam = 1:Nlam
    if (r==0)
        ialph=1; N_IL_MS(ilam,ialph) = 1; ialphst = 2;
    else
        ialphst = 1;
    end

    for ialph = ialphst:Nalph
        l = Lam(ilam,ialph); lam = l;
        a = Alph(ilam,ialph);

        filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
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
        0, 1, 1; ... % Green for 2 (ST)
        1, 0, 1; ... % Magenta for 3 (W)
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

figure;
x = [0 1]; y = [0 3];
imagesc(x,y,flipud(flagid));
set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
colormap(cMap);
cb = colorbar; caxis([0 4]);
xlabel('\alpha'); ylabel('\Lambda');
cb.Ticks = linspace(0,4,5); 
cb.TickLabels{1} = 'None';
cb.TickLabels{2} = 'IL';
cb.TickLabels{3} = 'ST';
cb.TickLabels{4} = 'W';
cb.TickLabels{5} = 'IL,W';
title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
saveas(gcf,[save_fig_dir figname '_imsc'],'png');

