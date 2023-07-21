%% This code plots stability boundary curves Lambda_i as a function of alph for different K and rho
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

Karr = [0 0.1 0.2 0.5 0.8 1 2]; NK = size(Karr,2);
rhoarr = [0.01 0.02 0.05 0.08 0.1]; Nrho = size(rhoarr,2);

% For reading files
Kstrarr = ["0" "0_1" "0_2" "0_5" "0_8" "1_0" "2_0"]; %"4_0"
rhostrarr = ["0_01" "0_02" "0_05" "0_08" "0_1"];

% For plotting
Kstrarrplot = ["0" "0.1" "0.2" "0.5" "0.8" "1" "2"]; %"4"
rhostrarrplot = ["0.01" "0.02" "0.05" "0.08" "0.1"];
colorarr = [0 0.4470 0.7410; 
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];
markarr = ["-" "--" "-." ":" "-" "--" "-."];
% markarr = ["-" "-" "-" "-" "-" "-" "-"]; % All solid lines

fdir1 = './files_Stab_boundary1/';
fdir2 = './files_Stab_boundary2/';
fdir3 = './files_Stab_boundary3/';

%% Plot Lam1,2,3 vs alph curve for different rho, fixed K

KstrarrB = ["0" "0_1" "0_2" "0_5" "0_8" "1" "2"]; %"4_0"

iK = 6; 
kappa = Karr(iK); Kstr = KstrarrB(iK);

% For saving figures
save_fig_dir = './figures_Lami/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname3 = sprintf('%s%s', 'Lami_v_alph_Td_K', Kstr);
figname4 = sprintf('%s%s', 'Lam12_v_alph_Td_K', Kstr);

idx = ones(size(rhoarr));
f3=figure; f4=figure;
for irho=1:Nrho
    r = rhoarr(irho) 
    rhostr = rhostrarr(irho);

    % Read files for Lam_1 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir1, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' rhostrarr(irho)]);
        idx(irho)=0;
        continue;
    end
    figure(f3); plot(alpharr,Lami,markarr(irho),'Color',colorarr(irho,:)); hold on;
    figure(f4); plot(alpharr,Lami,markarr(irho),'Color',colorarr(irho,:)); hold on;

end

% Lidx = logical(idx);
% Kstrarrplot2 = Kstrarrplot(Lidx);
figure(f3);  xlabel('\alpha'); ylabel('\Lambda_i');
% legend(Kstrarrplot2);


disp('Set 2 ********************************** ')
idx2 = ones(size(rhoarr));
for irho=1:Nrho
    r = rhoarr(irho) 
    rhostr = rhostrarr(irho);

    % Read files for Lam_2 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir2, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' rhostrarr(irho)]);
        idx2(irho)=0;
        continue;
    end
    figure(f3); plot(alpharr,Lami,markarr(irho),'Color',colorarr(irho,:)); hold on;
    figure(f4); plot(alpharr,Lami,markarr(irho),'Color',colorarr(irho,:),'HandleVisibility','off'); hold on;

end
figure(f4); xlabel('\alpha'); ylabel('\Lambda_1,\Lambda_2');
Lidx = logical(idx2);
rhostrarrplot2 = rhostrarrplot(Lidx);
legend(rhostrarrplot2);
disp('END OF Set 2 *************************** ');

idx3 = ones(size(rhoarr));
for irho=1:Nrho
    r = rhoarr(irho) 
    rhostr = rhostrarr(irho);

    % Read files for Lam_3 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir3, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' rhostrarr(irho)]);
        idx3(irho)=0;
        continue;
    end
    figure(f3); plot(alpharr,Lami,markarr(irho),'Color',colorarr(irho,:)); hold on;
end

figure(f3); xlim([0 1]); ylim([0 3]);
saveas(gcf,[save_fig_dir figname3],'png');

figure(f4); xlim([0 1]); ylim([0 1]);
saveas(gcf,[save_fig_dir figname4],'png');

disp('done - paused - exit!'); pause




%% Plot Lam1,2,3 vs alph curve for different K, fixed rho

irho = 3; 
r = rhoarr(irho); rhostr = rhostrarr(irho);

% For saving figures
save_fig_dir = './figures_Lami/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname1 = sprintf('%s%s', 'Lami_v_alph_Td_rho', rhostr);
figname2 = sprintf('%s%s', 'Lam2_v_alph_Td_rho', rhostr);

idx = ones(size(Karr));
f1=figure;
for iK=1:NK
    kappa = Karr(iK) 
    Kstr = Kstrarr(iK);

    % Read files for Lam_1 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir1, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' Kstrarr(iK)]);
        idx(iK)=0;
        continue;
    end
    figure(f1); plot(alpharr,smooth(Lami),markarr(iK),'Color',colorarr(iK,:)); hold on;

end

% Lidx = logical(idx);
% Kstrarrplot2 = Kstrarrplot(Lidx);
figure(f1);  xlabel('\alpha'); ylabel('\Lambda_i');
% legend(Kstrarrplot2);

disp('Set 2 ********************************** ')
idx2 = ones(size(Karr));
f2=figure;
for iK=1:NK
    kappa = Karr(iK) 
    Kstr = Kstrarr(iK);

    % Read files for Lam_2 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir2, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' Kstrarr(iK)]);
        idx2(iK)=0;
        continue;
    end
    figure(f1); plot(alpharr,smooth(Lami),markarr(iK),'Color',colorarr(iK,:)); hold on;
    figure(f2); plot(alpharr,smooth(Lami),markarr(iK),'Color',colorarr(iK,:)); hold on;

end
figure(f2); xlabel('\alpha'); ylabel('\Lambda_2');
Lidx = logical(idx2);
Kstrarrplot2 = Kstrarrplot(Lidx);
legend(Kstrarrplot2);
saveas(gcf,[save_fig_dir figname2],'png');
disp('END OF Set 2 *************************** ');


idx3 = ones(size(Karr));
for iK=1:NK
    kappa = Karr(iK) 
    Kstr = Kstrarr(iK);

    % Read files for Lam_3 boundary curve
    fname1 = sprintf('%s%s%s_%s%s', fdir3, 'Stab_bound_rho', rhostr, 'K', Kstr);
    filename = sprintf('%s%s',fname1,'.mat');
    if exist(filename,'file')
        load(filename);
    else
        disp(['----- Skipping' Kstrarr(iK)]);
        idx3(iK)=0;
        continue;
    end
    figure(f1); plot(alpharr,smooth(Lami),markarr(iK),'Color',colorarr(iK,:)); hold on;
end


figure(f1); xlim([0 1]); ylim([0 3]);
saveas(gcf,[save_fig_dir figname1],'png');

