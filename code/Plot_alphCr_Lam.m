%% This code plots the alph_cr as a function of Lambda
clear all; clc; close all
format short

FTsz = 26; 
% set(groot,'defaulttextFontName','Arial');
% set(groot,'defaultLegendFontName','Arial');
% set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

set(groot, 'DefaultLineLineWidth', 3.0);

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.1; rhostr = '0_1'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
ps=1; % self propulsion ratio: v01/v02
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

tol2=1e-2;
tol4 = 0.06; % tol for wall facing angle %0.14 for 1Td 0.2 case

%% Lambda
lamarr = [0:0.1:1.0 2.0 3.0 4.0]; Nlam = size(lamarr,2);
lamstrarr = ["0" "0_1" "0_2" "0_3" "0_4" "0_5" "0_6" "0_7" "0_8" "0_9" "1_0" ...
             "2_0" "3_0" "4_0"];
% lamstrarr2 = ["0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0" ...
%              "2.0" "3.0" "4.0"];

% lamarr = [0 0.1 0.2 0.5 1.0 2.0]; Nlam = size(lamarr,2);
% lamstrarr = ["0" "0_1" "0_2" "0_5" "1_0" "2_0"];
% lamstrarr2 = ["0" "0.1" "0.2" "0.5" "1.0" "2.0"];

width = 620; height = 400;
figure; set(gcf,'position',[161 415 width height]); 

Karr = [0.0:0.1:2.4 2.4+0.2:0.2:4]; NK = size(Karr,2); 

Karrplot = [0 0.1 0.2 0.5 0.8 1.0 2.0];
Kstrarr = ["0" "0_1" "0_2" "0_5" "0_8" "1_0" "2_0"];
Kstrarr2 = ["0" "0.1" "0.2" "0.5" "0.8" "1.0" "2.0"];

alpharr = 0.0:0.002:0.14; Nalph = size(alpharr,2);

colorarr = [0 0.4470 0.7410; 
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];
markarr = ["-" "--" "-." ":" "-" "--" "-."];

% Save figures
save_fig_dir = './figures_alphcr_Lam/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname1 = sprintf('%s%s', 'alphcr_v_K_finer_Tdrho', rhostr);    

for iK = 1:NK
    iK2= find(abs(Karr(iK)-Karrplot)<1e-5);
    if (isempty(iK2))
        continue;
    end
    kappa = Karr(iK);
    Kstr = Kstrarr(iK2);

    alphcr = nan(Nlam,1);

    for ilam = 1:Nlam
        l = lamarr(ilam);
        lam = l; Lamstr = lamstrarr(ilam); 

        % Read files
        fdir1 = './files_Roots_Kalph_finer/';
        fname = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'Lam', Lamstr);
       
        for ialph = 1:Nalph
            %%%%% Read roots from file ----------------------
            filename = sprintf('%s_%s%d_%s%d%s',fname,'iK', iK, 'ialph', ialph, '.mat');
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
    
                if (sum(idxS)+sum(idxAS)>0)
                    alphcr(ilam) = a;
                    break
                else
                    continue
                end
        
            end
    
        end
        disp(['Lam = ' num2str(l) ', alpha_cr = ' num2str(alphcr(ilam))]);
    end

    plot(alphcr,lamarr,markarr(iK2),'Color',colorarr(iK2,:)); hold on;
    xlabel('\alpha_{cr}'); ylabel('\Lambda');
    box on; %ylim([-0.5 0.5]);

    disp([' * * * * * * * Done with K = ' num2str(kappa) ' * * * * * * * ']);
end

legend(Kstrarr2);
xlim([0 0.2]); ylim([0 3]);
saveas(gcf,[save_fig_dir figname1],'png');
