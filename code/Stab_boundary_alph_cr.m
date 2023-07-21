%% This code finds the alph cr values for every lambda
clear all; clc; close all
format short

global r kappa a l  ps kappas kappap kappav

%% non dimensional parameters
r=0.01; rhostr = '0_01'; %dipole length rho
kappa=0; Kstr = '0'; %lateral line

ps=1; % self propulsion ratio: v01/v02
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

tol = 1e-5; 
tol2=1e-4; 

fdir1 = './files_Stab_bound_alphcr/';
if exist(fdir1, 'dir')==0
    mkdir(fdir1);
end

fname1 = sprintf('%s%s%s_%s%s', fdir1, 'Stab_bound_rho', rhostr, 'K', Kstr);

%% Solve - 5D eqns for 4 variables - specifying lambda
lamarr = [0.0:0.02:1.5 1.5+0.1:0.1:3.0]; Nlam = size(lamarr,2);
alpharr = 0.0:0.005:0.12; Nalph = size(alpharr,2);
[Alph,Lam] = meshgrid(alpharr,lamarr);

th10 = pi-pi/3;
th20 = [-pi/4 -2*pi/3]; nth20=size(th20,2);
th10B = [-pi/2 pi/2 -pi/2 pi/2];
th20B = [-pi/2+0.001 pi/2+0.001 pi/2+0.001 -pi/2+0.001]; nth0B=size(th20B,2);

% To write stab boundary curve in 
alphcr = nan(size(lamarr)); 

for ilam = 1:Nlam
    for ialph = 1:Nalph
        l = Lam(ilam,ialph); lam = l;
        a = Alph(ilam,ialph);

        ii = 0; clear RootsAll
        for zi10 = -0.5:0.02:0.5
            for zi20 = -0.5:0.02:0.5
                % initial conditions for numerical solver fsolve [zi1 zi2 th1 th2 lam] 
                %%%% zi20 sh be different from zi10 when lam = 0 to avoid NAN/INF
                zi20adj = zi20 + 0.001*(zi10==zi20);
                for ith20 = 1:nth20
                    X0 = [zi10 zi20adj th10 th20(ith20)];   
                    opts = optimset('Display','off');
                    opts.Algorithm = 'levenberg-marquardt';
                    [Root1,fval] = fsolve(@func_SOLVE5eq_var4_Tdip,X0,opts);
        %             zi1 = Root1(1); zi2 = Root1(2); th1 = Root1(3); th2 = Root1(4);
                    if (max(abs(fval))<tol) % Check that eqns are indeed satisfied
                        ii = ii + 1;
                        RootsAll(ii,:) = [Root1 lam]; % store all eq pts as row vectors
                    end
                end
            end
        end
        for zi10 = -0.475:0.475*2:0.475
            for zi20 = -0.475:0.475*2:0.475
                % initial conditions for numerical solver fsolve [zi1 zi2 th1 th2 lam] 
                %%%% zi20 sh be different from zi10 when lam = 0 to avoid NAN/INF
                zi20adj = zi20 + 0.001*(zi10==zi20);
                for ith0 = 1:nth0B
                    X0 = [zi10 zi20adj th10B(ith0) th20B(ith0)]; 
                    opts = optimset('Display','off');
                    opts.Algorithm = 'levenberg-marquardt';
                    [Root1,fval] = fsolve(@func_SOLVE5eq_var4_Tdip,X0,opts);
                    if (max(abs(fval))<tol) % Check that eqns are indeed satisfied
                        ii = ii + 1;
                        RootsAll(ii,:) = [Root1 lam]; % store all eq pts as row vectors
                    end
                end
            end
        end
        Roots = uniqueroots(RootsAll,tol);
        nuniq = size(Roots,1);
        eigflag = zeros(nuniq,1);
        E = zeros(nuniq,1);
        par = [r kappa  a ps kappap kappav  kappas];
        for irow=1:nuniq
            zi1 = Roots(irow,1); zi2 = Roots(irow,2); 
            th1 = Roots(irow,3); th2 = Roots(irow,4);
            lam = Roots(irow,5);
            x = [zi1 zi2 lam th1 th2];
            E(irow)=equationmap_Tdip(x,par); % Do this correction separately - cheaper
            eigflag(irow)=jac_Tdip(x,par); % eigmap = -1:unstable, 0:marg. stable, 1:stable 
        end
        Roots3 = [Roots eigflag];
        Roots3 = Roots3(E==1,:); % Remove roots that don't satisfy the eqns
        eigflag = Roots3(:,6);

        %%%%%%%% Check if 
        if (isempty(Roots3))
            disp('No roots!');
            continue;
        else
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

    if (~isnan(alphcr(ilam)))
        disp(['Done lambda = ' num2str(l) ', alph_cr = ' num2str(alphcr(ilam))]);
    else
        disp(['Could NOT find alpha_cr for Lambda = '  num2str(l)]);
    end

end

%% Save !!!!!
filename = sprintf('%s%s',fname1,'.mat');
save(filename,'alphcr','lamarr','r','kappa','ps','kappas','kappav','kappap');
