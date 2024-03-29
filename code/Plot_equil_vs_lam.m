%% This code plots the crosstream equilibria as a function of Lambda (for given alpha)
clear all; clc; close all
format short

FTsz = 20; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

set(groot, 'DefaultLineLineWidth', 2.5);

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.1; rhostr = '0_1'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
ps=1; % self propulsion ratio: v01/v02
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

% Pick which alpha you want
alphstr = '0';
a=0;

tol2=1e-2;
tol4 = 0.32; % tol for wall facing angle %0.14 for 1Td 0.22 case
% tol4B = 0.23; % tol for wall facing angle %0.14 for 1Td 0.22 case

% lamarr = 0.0:0.02:3.0; Nlam = size(lamarr,2);
% alpharr = [0.0:0.001:0.2 0.2+0.02:0.02:0.6]; Nalph = size(alpharr,2);
% [Alph,Lam] = meshgrid(alpharr,lamarr);
% l = Lam(ilam,1); lam = l
lamarr = [0.0:0.005:1.5 1.5+0.02:0.02:3.0]; Nlam = size(lamarr,2);

% Read files
fdir1 = './files1Td_v6_fine2/Roots_Lam/';
fname = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);

% Save figures
save_fig_dir = './figures_xiLam/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname1 = sprintf('%s%s_%s%s_%s%s', 'zi_v_Lam_Td_rho', rhostr, 'K', Kstr, 'alph', alphstr);
figname2 = sprintf('%s%s_%s%s_%s%s', 'th_v_Lam_Td_rho', rhostr, 'K', Kstr, 'alph', alphstr);

nS=0; nU=0;
% f1 = figure; 
f2a= figure; f2b= figure; f2c= figure; f2d= figure; 
f3a= figure; f3b= figure; 
f4a= figure; f4b= figure; 

% define arrays to store zi of stable and unstable eq points
maxEq = 8;
zi_Ups_IL_S = zeros(maxEq,Nlam); zi_Ups_IL_AS = zeros(maxEq,Nlam); zi_Ups_IL_U = zeros(maxEq,Nlam); 
zi_Ups_ST_S = zeros(maxEq,Nlam); zi_Ups_ST_AS = zeros(maxEq,Nlam); zi_Ups_ST_U = zeros(maxEq,Nlam); 
zi_Downs_IL_S = zeros(maxEq,Nlam); zi_Downs_IL_AS = zeros(maxEq,Nlam); zi_Downs_IL_U = zeros(maxEq,Nlam); 
zi_Downs_ST_S = zeros(maxEq,Nlam); zi_Downs_ST_AS = zeros(maxEq,Nlam); zi_Downs_ST_U = zeros(maxEq,Nlam); 
zi_Walls_IL_S = zeros(maxEq,Nlam); zi_Walls_IL_AS = zeros(maxEq,Nlam); zi_Walls_IL_U = zeros(maxEq,Nlam); 
zi_Walls_ST_S = zeros(maxEq,Nlam); zi_Walls_ST_AS = zeros(maxEq,Nlam); zi_Walls_ST_U = zeros(maxEq,Nlam); 

th_Ups_IL_S = nan(maxEq,Nlam); th_Ups_IL_AS = nan(maxEq,Nlam); th_Ups_IL_U = nan(maxEq,Nlam); 
th_Ups_ST_S = nan(maxEq,Nlam); th_Ups_ST_AS = nan(maxEq,Nlam); th_Ups_ST_U = nan(maxEq,Nlam); 
th_Downs_IL_S = nan(maxEq,Nlam); th_Downs_IL_AS = nan(maxEq,Nlam); th_Downs_IL_U = nan(maxEq,Nlam); 
th_Downs_ST_S = nan(maxEq,Nlam); th_Downs_ST_AS = nan(maxEq,Nlam); th_Downs_ST_U = nan(maxEq,Nlam); 
th_Walls_IL_S = nan(maxEq,Nlam); th_Walls_IL_AS = nan(maxEq,Nlam); th_Walls_IL_U = nan(maxEq,Nlam); 
th_Walls_ST_S = nan(maxEq,Nlam); th_Walls_ST_AS = nan(maxEq,Nlam); th_Walls_ST_U = nan(maxEq,Nlam); 

for ilam = 1:Nlam
    l = lamarr(ilam);
    disp(l);

    %%%%% Read roots from file ----------------------
    filename = sprintf('%s_%s%d_%s%s%s',fname,'ilam', ilam, 'alph', alphstr, '.mat');
%     filename = sprintf('%s_%s%d_%s%d%s',fname,'ilam', ilam, 'ialph', ialph, '.mat');
    load(filename);
%     Roots3

    %%%%% Pick U and S equilibrium pts
    if (isempty(Roots3))
        disp('No roots!');
        continue;
    else
        eigflag = Roots3(:,6);

        %%%% eigflag = -1:U, 0:MS, 1:S 
        % Asymp stable (AS) 
        idxAS = eigflag > 0; 
        % Marg stable (S)
        idxS = eigflag == 0; %N_S(ilam,ialph) = sum(idxS);
        % Unstable (U)
        idxU = eigflag < 0; %N_U(ilam,ialph) = sum(idxU);
        % IL position
        idxIL = abs(Roots3(:,1)-Roots3(:,2))<tol2;
        % ST position (staggered symm about centerline)
        idxST = ( abs(Roots3(:,1)+Roots3(:,2))<tol2 ) & ...
                  abs(Roots3(:,1))>tol2 & abs(Roots3(:,2))>tol2;

        RootsIL_S = Roots3((idxS==1 & idxIL==1),:);
        RootsIL_AS = Roots3((idxAS==1 & idxIL==1),:); 
        RootsIL_U = Roots3((idxU==1 & idxIL==1),:);

        RootsST_S = Roots3((idxS==1 & idxST==1),:); 
        RootsST_AS = Roots3((idxAS==1 & idxST==1),:); 
        RootsST_U = Roots3((idxU==1 & idxST==1),:);

        % Swimming upstream - plot zi1
        idxIL_S = abs(abs(RootsIL_S(:,3))-pi)<tol4; 
        idxIL_AS = abs(abs(RootsIL_AS(:,3))-pi)<tol4; 
        idxIL_U = abs(abs(RootsIL_U(:,3))-pi)<tol4;

        idxST_S = abs(abs(RootsST_S(:,3))-pi)<tol4; 
        idxST_AS = abs(abs(RootsST_AS(:,3))-pi)<tol4; 
        idxST_U = abs(abs(RootsST_U(:,3))-pi)<tol4;

        count = sum(idxIL_S) + sum(idxIL_U) + sum(idxIL_AS) + sum(idxST_S) + sum(idxST_U) + sum(idxST_AS); 
        
        figure(f2a);
        plot(l*ones(size(RootsIL_S(idxIL_S,1))),RootsIL_S(idxIL_S,1),'xg'); hold on
        plot(l*ones(size(RootsIL_AS(idxIL_AS,1))),RootsIL_AS(idxIL_AS,1),'s','MarkerFaceColor',[0 0.5 0]); hold on
        plot(l*ones(size(RootsIL_U(idxIL_U,1))),RootsIL_U(idxIL_U,1),'or'); hold on
        title('Inline upstream');
        figure(f2b);
        plot(l*ones(size(RootsST_S(idxST_S,1))),RootsST_S(idxST_S,1),'og'); hold on
        plot(l*ones(size(RootsST_AS(idxST_AS,1))),RootsST_AS(idxST_AS,1),'s','MarkerFaceColor',[0 0.5 0]); hold on
        plot(l*ones(size(RootsST_U(idxST_U,1))),RootsST_U(idxST_U,1),'xr'); hold on
        title('Staggered upstream');

        if (sum(idxIL_S)>0)
            zi_Ups_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,1);
            th_Ups_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,3);
        end
        if (sum(idxIL_AS)>0)
            zi_Ups_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,1);
            th_Ups_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,3);
        end
        if (sum(idxIL_U)>0)
            zi_Ups_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,1);
            th_Ups_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,3);
        end
        if (sum(idxST_S)>0)
            zi_Ups_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,1);
            th_Ups_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,3);
        end
        if (sum(idxST_AS)>0)
            zi_Ups_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,1);
            th_Ups_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,3);
        end
        if (sum(idxST_U)>0)
            zi_Ups_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,1);
            th_Ups_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,3);
        end

        % Swimming downstream - plot zi1
        idxIL_S = abs(RootsIL_S(:,3))<tol4; 
        idxIL_AS = abs(RootsIL_AS(:,3))<tol4; 
        idxIL_U = abs(RootsIL_U(:,3))<tol4;
        
        idxST_S = abs(RootsST_S(:,3))<tol4; 
        idxST_AS = abs(RootsST_AS(:,3))<tol4;
        idxST_U = abs(RootsST_U(:,3))<tol4;
        
        count = count + sum(idxIL_S) + sum(idxIL_AS) + sum(idxIL_U) + sum(idxST_S) + sum(idxST_AS) + sum(idxST_U);
        
        figure(f3a);
        plot(l*ones(size(RootsIL_S(idxIL_S,1))),RootsIL_S(idxIL_S,1),'xg'); hold on
        plot(l*ones(size(RootsIL_AS(idxIL_AS,1))),RootsIL_AS(idxIL_AS,1),'s','Color',[0 0.5 0]); hold on
        plot(l*ones(size(RootsIL_U(idxIL_U,1))),RootsIL_U(idxIL_U,1),'xr'); hold on
        title('Inline downstream');
        figure(f3b);
        plot(l*ones(size(RootsST_S(idxST_S,1))),RootsST_S(idxST_S,1),'og'); hold on
        plot(l*ones(size(RootsST_AS(idxST_AS,1))),RootsST_AS(idxST_AS,1),'s','Color',[0 0.5 0]); hold on
        plot(l*ones(size(RootsST_U(idxST_U,1))),RootsST_U(idxST_U,1),'or'); 
        title('Staggered downstream');
        if (sum(idxIL_S)>0)
            zi_Downs_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,1);
            th_Downs_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,3);
        end
        if (sum(idxIL_AS)>0)
            zi_Downs_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,1);
            th_Downs_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,3);
        end
        if (sum(idxIL_U)>0)
            zi_Downs_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,1);
            th_Downs_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,3);
        end
        if (sum(idxST_S)>0)
            zi_Downs_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,1);
            th_Downs_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,3);
        end
        if (sum(idxST_AS)>0)
            zi_Downs_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,1);
            th_Downs_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,3);
        end
        if (sum(idxST_U)>0)
            zi_Downs_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,1);
            th_Downs_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,3);
        end
 
        % Swimming wall - plot zi1
        idxIL_S = abs(abs(RootsIL_S(:,3))-pi/2)<tol4; 
        idxIL_AS = abs(abs(RootsIL_AS(:,3))-pi/2)<tol4; 
        idxIL_U = abs(abs(RootsIL_U(:,3))-pi/2)<tol4;

        idxST_S = abs(abs(RootsST_S(:,3))-pi/2)<tol4; 
        idxST_AS = abs(abs(RootsST_AS(:,3))-pi/2)<tol4; 
        idxST_U = abs(abs(RootsST_U(:,3))-pi/2)<tol4;

        % Only consider points at walls
        idxIL_S_W = abs(RootsIL_S(:,1))>0.4; 
        idxIL_AS_W = abs(RootsIL_AS(:,1))>0.4; 
        idxIL_U_W = abs(RootsIL_U(:,1))>0.4; 

        idxST_S_W = abs(RootsST_S(:,1))>0.4; 
        idxST_AS_W = abs(RootsST_AS(:,1))>0.4; 
        idxST_U_W = abs(RootsST_U(:,1))>0.4; 

        idxIL_S = idxIL_S==1 & idxIL_S_W==1; 
        idxIL_AS = idxIL_AS==1 & idxIL_AS_W==1; 
        idxIL_U = idxIL_U==1 & idxIL_U_W==1;

        idxST_S = idxST_S==1 & idxST_S_W==1; 
        idxST_AS = idxST_AS==1 & idxST_AS_W==1;
        idxST_U = idxST_U==1 & idxST_U_W==1;

        count = count + sum(idxIL_S) + sum(idxIL_AS) + sum(idxIL_U) + ...
                        sum(idxST_S) + sum(idxST_AS) + sum(idxST_U);
        figure(f4a);
        plot(l*ones(size(RootsIL_S(idxIL_S,1))),RootsIL_S(idxIL_S,1),'xg'); hold on
        plot(l*ones(size(RootsIL_AS(idxIL_AS,1))),RootsIL_AS(idxIL_AS,1),'s','Color',[0 0.5 0]); hold on
        plot(l*ones(size(RootsIL_U(idxIL_U,1))),RootsIL_U(idxIL_U,1),'xr'); hold on
        title('Inline wall');
        figure(f4b);
        plot(l*ones(size(RootsST_S(idxST_S,1))),RootsST_S(idxST_S,1),'og'); hold on
        plot(l*ones(size(RootsST_AS(idxST_AS,1))),RootsST_AS(idxST_AS,1),'s','Color',[0 0.5 0]); hold on
        plot(l*ones(size(RootsST_U(idxST_U,1))),RootsST_U(idxST_U,1),'or'); 
        title('Staggered wall');
        if (sum(idxIL_S)>0)
            zi_Walls_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,1);
            th_Walls_IL_S(1:sum(idxIL_S),ilam) = RootsIL_S(idxIL_S,3);
        end
        if (sum(idxIL_AS)>0)
            zi_Walls_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,1);
            th_Walls_IL_AS(1:sum(idxIL_AS),ilam) = RootsIL_AS(idxIL_AS,3);
        end
        if (sum(idxIL_U)>0)
            zi_Walls_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,1);
            th_Walls_IL_U(1:sum(idxIL_U),ilam) = RootsIL_U(idxIL_U,3);
        end
        if (sum(idxST_S)>0)
            zi_Walls_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,1);
            th_Walls_ST_S(1:sum(idxST_S),ilam) = RootsST_S(idxST_S,3);
        end
        if (sum(idxST_AS)>0)
            zi_Walls_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,1);
            th_Walls_ST_AS(1:sum(idxST_AS),ilam) = RootsST_AS(idxST_AS,3);
        end
        if (sum(idxST_U)>0)
            zi_Walls_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,1);
            th_Walls_ST_U(1:sum(idxST_U),ilam) = RootsST_U(idxST_U,3);
        end

        if (count/size(Roots3,1)<1)
            disp('Error! Some config is not one of IL-ST-Wall');
            disp(size(Roots3,1)-count);%pause;
        end

    end

%     pause;

end

disp(['IL Ups S --> ' num2str(min(abs(th_Ups_IL_S(:)))/pi) 'pi : ' num2str(max(abs(th_Ups_IL_S(:)))/pi) 'pi' ...
      ',  U --> ' num2str(min(abs(th_Ups_IL_U(:)))/pi) 'pi : ' num2str(max(abs(th_Ups_IL_U(:)))/pi) 'pi']);

disp(['IL Downs S --> ' num2str(min(abs(th_Downs_IL_S(:)))/pi) 'pi' ':' num2str(max(abs(th_Downs_IL_S(:)))/pi)  'pi'...
      ',  U --> ' num2str(min(abs(th_Downs_IL_U(:)))/pi) 'pi' ':' num2str(max(abs(th_Downs_IL_U(:)))/pi) 'pi']);

disp(['ST Ups S --> ' num2str(min(abs(th_Ups_ST_S(:)))/pi) 'pi' ':' num2str(max(abs(th_Ups_ST_S(:)))/pi) 'pi' ...
      ',  U --> ' num2str(min(abs(th_Ups_ST_U(:)))/pi) 'pi' ':' num2str(max(abs(th_Ups_ST_U(:)))/pi) 'pi']);

disp(['ST Downs S --> ' num2str(min(abs(th_Downs_ST_S(:)))/pi) 'pi' ':' num2str(max(abs(th_Downs_ST_S(:)))/pi) 'pi' ...
      ',  U --> ' num2str(min(abs(th_Downs_ST_U(:)))/pi) 'pi' ':' num2str(max(abs(th_Downs_ST_U(:)))/pi) 'pi']);

disp(['IL Walls AS --> ' num2str(min(abs(th_Walls_IL_AS(:)))/pi) 'pi' ':' num2str(max(abs(th_Walls_IL_AS(:)))/pi) 'pi' ...
      ',  U --> ' num2str(min(abs(th_Walls_IL_U(:)))/pi) 'pi' ':' num2str(max(abs(th_Walls_IL_U(:)))/pi) 'pi']);

disp(['ST Walls AS --> ' num2str(min(abs(th_Walls_ST_AS(:)))/pi) 'pi' ':' num2str(max(abs(th_Walls_ST_AS(:)))/pi) 'pi' ...
      ',  U --> ' num2str(min(abs(th_Walls_ST_U(:)))/pi) 'pi' ':' num2str(max(abs(th_Walls_ST_U(:)))/pi) 'pi']);


% height= 1293; width = 420;
width = 420; height = 200;

% figure(f1); set(gcf,'position',[961 415 width height]); box on; ylim([-0.5 0.5]);
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_All'],'png');
% 
figure(f2a); set(gcf,'position',[400 830 width height]); box on;  ylim([-0.5 0.5]);
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_IL_Ups'],'png');
% 
figure(f2b); set(gcf,'position',[400 830 width height]); box on; ylim([-0.5 0.5]); 
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_ST_Ups'],'png');

figure(f3a); set(gcf,'position',[400 830 width height]); box on; ylim([-0.5 0.5]); 
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_IL_Downs'],'png');
% 
figure(f3b); set(gcf,'position',[400 830 width height]); box on; ylim([-0.5 0.5]); 
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_ST_Downs'],'png');
% 
figure(f4a); set(gcf,'position',[400 830 width height]); box on; ylim([-0.5 0.5]); 
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_IL_Wall'],'png');
% 
figure(f4b); set(gcf,'position',[400 830 width height]); box on; ylim([-0.5 0.5]); 
% xlabel('\alpha'); ylabel('\xi'); saveas(gcf,[save_fig_dir figname1 '_ST_Wall'],'png');

figure; set(gcf,'position',[161 415 width height]);
plot_lines2(zi_Ups_IL_U,lamarr,-1,0.06);
plot_lines2(zi_Ups_IL_S,lamarr,0,0.06);
plot_lines2(zi_Ups_IL_AS,lamarr,1,0.06);
% plot(lamarr(1:15),zeros(size(lamarr(1:15))),'-r');
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_IL_Ups'],'png');

figure; set(gcf,'position',[561 415 width height]);
plot_lines2(zi_Ups_ST_U,lamarr,-1,0.05);
plot_lines2(zi_Ups_ST_S,lamarr,0,0.05);
plot_lines2(zi_Ups_ST_AS,lamarr,1,0.05);
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_ST_Ups'],'png');

figure; set(gcf,'position',[961 415 width height]);
plot_lines2(zi_Downs_IL_U,lamarr,-1,0.08);
plot_lines2(zi_Downs_IL_S,lamarr,0,0.08);
plot_lines2(zi_Downs_IL_AS,lamarr,1,0.08);
% plot(lamarr(1:15),zeros(size(lamarr(1:15))),'-r');
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_IL_Downs'],'png');

figure; set(gcf,'position',[161 115 width height]);
plot_lines2(zi_Downs_ST_U,lamarr,-1,0.08);
plot_lines2(zi_Downs_ST_S,lamarr,0,0.08);
plot_lines2(zi_Downs_ST_AS,lamarr,0,0.08);
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_ST_Downs'],'png');

figure; set(gcf,'position',[561 115 width height]);
plot_lines2(zi_Walls_IL_U,lamarr,-1,0.08);
plot_lines2(zi_Walls_IL_S,lamarr,0,0.08);
plot_lines2(zi_Walls_IL_AS,lamarr,1,0.08);
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_IL_Walls'],'png');

figure; set(gcf,'position',[961 115 width height]);
plot_lines2(zi_Walls_ST_U,lamarr,-1,0.08);
plot_lines2(zi_Walls_ST_S,lamarr,0,0.08);
plot_lines2(zi_Walls_ST_AS,lamarr,1,0.08);
box on; ylim([-0.5 0.5]); xlim([0 3.0]); xlabel('\Lambda'); ylabel('\xi'); 
saveas(gcf,[save_fig_dir figname1 '_line_ST_Walls'],'png');


% disp('ADD A RED LINE for Ups_IL_U & Downs_IL_U at centerline');
% disp('those two are equilibria at the missing low lambda values, I checked!');
% disp('Add AS - asymp stable - dark green color');