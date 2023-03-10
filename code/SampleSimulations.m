% Sample simulations: Numerical integration of the 5D governing equations 
% of the fish dynamics to show time-evolution of the 
% dimensionless position variables & orientations of both fish

clear all,close all, clc  
FTsz = 20; 
set(groot,'defaulttextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);
%global variables linked to the function ODEfive
global alpha k rho  psi kappap kappav kappas

%non dimensional parameters
k=0; % lateral line gain
rho=0.1; %dipole length
alpha=0.4; %flow speed/fish speed ratio
psi=1; % self propulsion ratio v02/v01
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

%% Initial conditions: 

lambda0 = 1.0; %inter-fish distance (stream-wise)
zi10=0.00+0.02*randn(1); %cross-stream coordinates
zi20=-0.00+0.02*randn(1);
thetaf10=pi+0.04*randn(1); %heading anglese
thetaf20=pi+0.04*randn(1);

y0 = [zi10 zi20 lambda0 thetaf10 thetaf20]'; %(5x1 vector with intial contidions)

%% ODE INTEGRATION using ode113 (better for non-linear systems)

tf=100; dt = 0.001; % time of integration (can be changed)
tspan = [0 tf]; 
tstart = tic;
tarr = 0:dt:tf;
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'Stats','on','OutputFcn',@odeplot,...
            'Events',@(t,y) myevent(t,y,tstart));%set tolerances
[t,y] = ode113(@ODEfive_Tdip,tarr,y0); %integration of the equations contained in ODEfive.m function
   
iRoots = abs(imag(y));tol=1e-4;
if (max(max(iRoots))>tol)
    for it=1:size(t,1)
        if (max(abs(imag(y(it,:))))>tol)
            itstop=it-1;
            disp('imag!'); y0'
            break;
        end
    end
    t = t(1:itstop); 
    y = y(1:itstop,:); y = real(y);
else
    y = real(y);
end

%% Plot
figure;
plot(t,y(:,1),'r',t,y(:,2),'b','LineWidth',1.5);
xlabel('Time (s)','FontSize', 20),ylabel('\xi','FontSize', 20);
legend('\xi_1','\xi_2');
ylim([-0.5 0.5]);pbaspect([1.65 1 1])

figure;
plot(t,y(:,4),'r',t,y(:,5),'b','LineWidth',1.5);
xlabel('Time (s)','FontSize', 20),ylabel('\theta','FontSize', 20);
legend('\theta_1','\theta_2');
ylim([pi/2 3*pi/2]); 
yticks([pi/2 pi 3*pi/2]); yticklabels({'\pi/2','\pi','3\pi/2'});
pbaspect([1.65 1 1]);

figure; %stream-wise
plot(t,y(:,3),'b','LineWidth',1.5)
xlabel('Time (s)','FontSize', 20),ylabel('\Lambda','FontSize', 20)
pbaspect([1.65 1 1])



%% Function
function [values,isterminal,direction] = myevent(t,y,tstart)
 %  Don't let t cross zero...a dummy "event" to illustrate how 
 %  one might handle other events in conjunction with the time
 %  constraint.  Not necessary, but I put it in just in case.
 values(1) = t;
 %  Don't let integration go for more than 1.2 seconds.
 values(2) = toc(tstart) < 60;
 isterminal = true(size(values));
 direction = zeros(size(values));
end