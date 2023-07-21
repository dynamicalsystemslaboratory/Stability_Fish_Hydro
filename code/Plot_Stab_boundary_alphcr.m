%% This code plots the alph_cr as a function of Lambda 
% for different K (fixed rho)
% and for different rho (fixed K)

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
Kstrarr = ["0" "0_1" "0_2" "0_5" "0_8" "1" "2"]; %"4_0"
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

fdir1 = './files_Stab_bound_alphcr/';

%% Plot Lam1,2,3 vs alph curve for different rho, fixed K

KstrarrB = ["0" "0_1" "0_2" "0_5" "0_8" "1" "2"]; %"4_0"

iK = 1; 
kappa = Karr(iK); Kstr = KstrarrB(iK);

% For saving figures
save_fig_dir = './figures_alphcr/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname3 = sprintf('%s%s', 'Lam_v_alphcr_Td_K', Kstr);

idx = ones(size(rhoarr));
f3=figure; 
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
%     pause

    % Make sure we only take alphcr curve
    Nlam = size(lamarr,2);
    for ilam=Nlam:-1:1
        if (alphcr(ilam)==0 && max(abs(alphcr))>0 ) % && lamarr(ilam)<1.0)
            alphcr(ilam)=nan;
        end
    end

%     alphcr2 = alphcr;
%     alphcr2 = [smooth(alphcr(1:75))' alphcr(76:end)];

    alphcr2 = [alphcr(1:4:60) alphcr(62:2:end)];
    lamarr2 = [lamarr(1:4:60) lamarr(62:2:end)];
%     alphcr2 = alphcr(1:3:end);
%     lamarr2 = lamarr(1:3:end);
    figure(f3); plot(alphcr2,lamarr2,markarr(irho),'Color',colorarr(irho,:)); hold on;
end

figure(f3);  xlabel('\alpha_{cr}'); ylabel('\Lambda');
Lidx = logical(idx);
rhostrarrplot2 = rhostrarrplot(Lidx);
legend(rhostrarrplot2);
xlim([0 0.2]); ylim([0 3]);
saveas(gcf,[save_fig_dir figname3],'png');

disp('done - paused - exit!'); pause

%% Plot Lam1,2,3 vs alph curve for different K, fixed rho

irho = 5; 
r = rhoarr(irho); rhostr = rhostrarr(irho);

% For saving figures
save_fig_dir = './figures_alphcr/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end
figname1 = sprintf('%s%s', 'Lam_v_alphcr_Td_rho', rhostr);

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

    % Make sure we only take alphcr curve
    Nlam = size(lamarr,2);
    for ilam=Nlam:-1:1
        if (alphcr(ilam)==0 && max(abs(alphcr))>0 ) % && lamarr(ilam)<1.0)
            alphcr(ilam)=nan;
        end
    end
%     alphcr2 = alphcr;
%     alphcr2 = [smooth(alphcr(1:75))' alphcr(76:end)];
    alphcr2 = [alphcr(1:4:60) alphcr(62:2:end)];
    lamarr2 = [lamarr(1:4:60) lamarr(62:2:end)];
    figure(f1); plot(alphcr2,lamarr2,markarr(iK),'Color',colorarr(iK,:)); hold on;

%     pause
end

Lidx = logical(idx);
Kstrarrplot2 = Kstrarrplot(Lidx);
figure(f1);  xlabel('\alpha_{cr}'); ylabel('\Lambda');
legend(Kstrarrplot2);
xlim([0 0.2]); ylim([0 3]);
saveas(gcf,[save_fig_dir figname1],'png');

