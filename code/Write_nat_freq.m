%% This code computes the natural frequencies of stable equilibria
% at different values of alpha-lambda
clear all; clc; close all
format short

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.1; rhostr = '0_1'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line
tol = 1e-05; tol2=1e-3; 
tol4 = 0.1; % angle tol for IL
tol5 = 0.8; % angle tol for ST 
tol6 = 0.2; % angle tol for WF 
% When doing for rho=0 case, remove ialph=1 values!

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

om1 = nan(Nlam,Nalph);om2 = nan(Nlam,Nalph);

% Folder Roots
lamarr_coarse = 0.0:0.02:3.0; Nlam_coarse = size(lamarr_coarse,2);
alpharr_coarse = 0.0:0.02:1.0; Nalph_coarse = size(alpharr_coarse,2);
Nlam1_coarse = size(0:0.02:1.5,2); Nalph1_coarse = size(0:0.02:0.16,2);

fdir = './files1Td_v6_fine2/';
fdirB = './../MATLAB_CODES/files1Td_v6_fine2/';

% Read in parameters from one of the roots files
ilam = 1; ialph = 1; 
fdir1 = [fdir 'Roots3/'];
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);
filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
load(filename);

fdir2 = [fdir 'om_fine/'];
if exist(fdir2, 'dir')==0
    mkdir(fdir2);
end
fname = sprintf('%s%s_%s_%s_%s', fdir2, 'w1w2_rho', rhostr, 'K', Kstr);


%%
for ilam = 1:Nlam
    for ialph = 2:Nalph
        disp(ialph)
        l = Lam(ilam,ialph); lam = l;
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


        if (isempty(Roots3))
            continue;
        else
            eigflag = Roots3(:,6);

            %%%% eigflag = -1:U, 0:MS, 1:S 
            % Stable (S) or Marg stable (MS)
            idxS = eigflag >= 0; 
            Roots = Roots3(idxS,:); nrow = size(Roots,1);
            
            if (isempty(Roots))
                continue;
            else
                om_all=[];
                for irow=1:nrow
                    zi1 = Roots(irow,1); zi2 = Roots(irow,2); 
                    th1 = Roots(irow,3); th2 = Roots(irow,4);
                    lam = Roots(irow,5);
                    x = [zi1 zi2 lam th1 th2];
                    par = [r kappa  a ps kappap kappav  kappas];
        %             equationmap_Tdip(x,par)
                    Jvectmp = func_jacob_Tdip(x,par); % store all Jacobians at these eq pts as row vectors
                    J = reshape(Jvectmp,[5,5]);  %% Jac order is [zi1 zi2 th1 th2 lam] 
                    [V,D]=eig(J);
                    om = abs(imag(diag(D)));
                    om = uniquetol(om, 0.01, 'DataScale', 1,'Byrows',true);
                    om = om(om>0.001); % keep only non zero ones
                    if (size(om)>2)
                        disp('===================== ERROR !!! ==========');
                        disp(D);
                    end
                    om_all = [om_all; om];    
                end
                om_all = uniquetol(om_all, 0.01, 'DataScale', 1,'Byrows',true);

                if (size(om_all)>2)
                    disp(['Error! Multiple S pts - different om1,2 - at lam = ' num2str(l) ', alph = ' num2str(a)]); 
                    disp(om_all);
                elseif (~isempty(om_all))
                    om1(ilam,ialph) = min(om_all); om2(ilam,ialph) = max(om_all);
                end

            end

        end

    end
    disp(['Done lambda = ' num2str(l)]);
end

filename = sprintf('%s%s',fname,'_B.mat');
save(filename,'om1','om2','Lam', 'Alph','r','kappa','ps','kappas','kappav','kappap');


