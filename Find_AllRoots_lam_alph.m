%% Find the unique roots of the system of equations

clear all; close all
format short
global r kappa a l  ps kappas kappap kappav
%% non dimensional parameters
r=0.1; rhostr = '0_1'; %dipole length rho
% a=0.15; %flow speed/fish speed ratio
ps=1; % self propulsion ratio: v01/v02
kappa=0; Kstr = '0'; %lateral line

fdir = './files1Td_v6_fine2/';
if exist(fdir, 'dir')==0
    mkdir(fdir);
end

fdir1 = './files1Td_v6_fine2/Roots/';
if exist(fdir1, 'dir')==0
    mkdir(fdir1);
end

fname = sprintf('%s%s_%s_%s_%s%s', fdir, 'Num_rho', rhostr, 'K', Kstr, '.mat');
fname1 = sprintf('%s%s_%s_%s_%s', fdir1, 'Roots_rho', rhostr, 'K', Kstr);

% Vision
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

tol = 1e-5; 
tol2 = 1e-4; % For checking if near enough to be tandem, etc.

%% Solve - 5D eqns for 4 variables - specifying lambda
lamarr = 0.0:0.02:3.0; Nlam = size(lamarr,2);
alpharr = 0.0:0.02:1.0; Nalph = size(alpharr,2);
% lamarr = 0.0:0.5:1.0; Nlam = size(lamarr,2);
% alpharr = 0.0:0.5:1.0; Nalph = size(alpharr,2);
[Alph,Lam] = meshgrid(alpharr,lamarr);

th10 = -3*pi/4:pi/3:11*pi/12; th10 = [th10 pi/2 -pi/2];
th20 = -3*pi/4:pi/3:11*pi/12; th20 = [th20 pi/2+0.001 -pi/2+0.001];
nth10=size(th10,2); nth20=size(th20,2);

for ilam = 1:10
    for ialph = 1:Nalph
%         disp(ialph)
        l = Lam(ilam,ialph); lam = l;
        a = Alph(ilam,ialph);

        ii = 0; clear RootsAll
        for zi10 = -0.5:0.02:0.5
            for zi20 = -0.5:0.02:0.5
                % initial conditions for numerical solver fsolve [zi1 zi2 th1 th2 lam] 
                %%%% zi20 sh be different from zi10 when lam = 0 to avoid NAN/INF
                zi20adj = zi20 + 0.001*(zi10==zi20);
                for ith10 = 1:nth10
                    for ith20 = 1:nth20
                        X0 = [zi10 zi20adj th10(ith10) th20(ith20)];   
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

        % Save all the roots for given rho val K val (so we can calc Jac from roots later)!!!!!
        filename = sprintf('%s_%s%d_%s%d%s',fname1,'ilam', ilam, 'ialph', ialph, '.mat');
        save(filename,'Roots3','l','a','r','kappa','ps','kappas','kappav','kappap');

    end
    disp(['Done lambda = ' num2str(l)]);

end

