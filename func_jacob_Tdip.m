function Avec = func_jacob_Tdip(y,par)
% This function returns the 5x5 Jacobian matrix
% in a vectorized form

% global  r kappa  a ps kappap kappav kappas
r = par(1); kappa = par(2); a = par(3); ps = par(4); 
kappap = par(5); kappav = par(6); kappas = par(7);

syms zita1 zita2 thetaf1 thetaf2 rho h  k lambda L
psi =sym("psi");
ks=sym("ks");
kp=sym("kp");
kv=sym("kv");
alpha=sym("alpha");

load('Jac_Tdip.mat','A');

A1 = subs(A,[rho,k,alpha,psi,ks,kp,kv], [r,kappa,a,ps,kappas,kappap,kappav]); %substitute parameters
A2 = double(subs(A1,[zita1,zita2,lambda,thetaf1,thetaf2], [y(1),y(2),y(3),y(4),y(5)])); %substitute variables

Avec = reshape(A2,[1,5*5]);

end