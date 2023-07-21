%% This code writes the component of zero eigenvector of stable equilibria along lambda
clear all

FTsz = 20; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

%% Roots3
global r kappa a l  ps kappas kappap kappav

% ilam = 18; ialph=45;  
r=0.1; rhostr = '0_1'; %dipole length rho
% a=0.15; %flow speed/fish speed ratio
ps=1; % self propulsion ratio: v01/v02
kappa=0; Kstr = '0'; %lateral line


fdir = './files1Td_v6_fine2/';
fdir1 = [fdir 'Roots/'];
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);

filename2 = sprintf('%s%s',fdir, 'Ezero_fine.mat'); % To save

save_fig_dir = './figures_natfreq/';

% fdir = './files4Td_v4_fine2/';
% fname2 = sprintf('%s%s_%s_%s_%s%s', fdir, 'Num_rho', rhostr, 'K', Kstr, '.mat');
% load(fname2);

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

fdir = './files1Td_v6_fine2/';
fdirB = './../MATLAB_CODES/files1Td_v6_fine2/';

tol=1e-6;

Ezero = nan(Nlam,Nalph);

for ilam = 1:Nlam
    l = Lam(ilam,1)
    for ialph=1:Nalph
        a = Alph(ilam,ialph);

%         filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
%         load(filename);
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

        idxS = Roots3(:,6)>=0; RootsS = Roots3(idxS,:); 

        if (~isempty(RootsS))
            for irow=1:size(Roots3,1)
                zi1 = Roots3(irow,1); zi2 = Roots3(irow,2); 
                th1 = Roots3(irow,3); th2 = Roots3(irow,4);
                lam = Roots3(irow,5);
                x = [zi1 zi2 lam th1 th2];
                par = [r kappa  a ps kappap kappav  kappas];
                Jvectmp = func_jacob_Tdip(x,par); % store all Jacobians at these eq pts as row vectors
                J = reshape(Jvectmp,[5,5]);  %% Jac order is [zi1 zi2 th1 th2 lam] 
                [V,D]=eig(J);
                d = diag(D);
                [dzero,i]= min(abs(d)); % zero eigenvalue
                vzero = V(:,i); % zero eigenvector
                idx = abs(imag(vzero))<tol;
                vzero(idx) = real(vzero(idx));
                Vzero(irow,:) = vzero';
            end
            % Stable
            VzeroS = Vzero(idxS,:); NS = size(RootsS,1);
            if (NS==1)
                Ezero(ilam,ialph)=abs(VzeroS(1,5));
            else
                Ezerotmp = uniquetol(abs(VzeroS(:,5)),0.01, 'DataScale', 1,'Byrows',true);
                if size(Ezerotmp,1)>1
                    disp('Error!'); disp(VzeroS); pause;
                else
                    Ezero(ilam,ialph)=Ezerotmp;
                end
            end

        end

    end
end

save(filename2,'Ezero','Lam','Alph');
% 
% % Plot
% 
% figure;
% [c,h] = contourf(Alph,Lam,Ezero); 
% colorbar;
% set(h, 'edgecolor','none');
% xlabel('\alpha'); ylabel('\Lambda');
% saveas(gcf,[save_fig_dir 'zeroEvec_lam'],'png');
% 
% figure;
% x = [0 1]; y = [0 3];
% imagesc(x,y,flipud(Ezero));
% set(gca,'YTick',get(gca, 'YTick'),'YTickLabel',[3 2.5 2 1.5 1 0.5 0]);
% % colormap(cMap);
% colorbar; 
% % caxis([0 4]);
% xlabel('\alpha'); ylabel('\Lambda');
% % title(['Stable: K = ' Kstr ' &  \rho = ' num2str(r) ' (h = ' num2str(0.2/r) ' BL)' ]);
% saveas(gcf,[save_fig_dir 'zeroEvec_lam_imsc'],'png');

