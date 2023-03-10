function EIG = jac_Tdip(y,par)

% This function classifies the equilibrium point as Stable, Unstable, and
% Marginally stable based on eigendecomposition of the Jacobian matrix

r = par(1); kappa = par(2); a = par(3); ps = par(4); 
kappap = par(5); kappav = par(6); kappas = par(7);

syms zita1 zita2 thetaf1 thetaf2 rho h  k lambda L
psi =sym("psi");
ks=sym("ks");
kp=sym("kp");
kv=sym("kv");
alpha=sym("alpha");

% Jacobian 5x5: state matrix
load('Jac_Tdip.mat','A');

A1 = subs(A,[rho,k,alpha,psi,ks,kp,kv], [r,kappa,a,ps,kappas,kappap,kappav]); %substitute parameters
A2 = double(subs(A1,[zita1,zita2,lambda,thetaf1,thetaf2], [y(1),y(2),y(3),y(4),y(5)])); %substitute variables
% [V,D]=eig(A2);
eigenval = eigs(A2,5,'largestreal'); % sort from highest to lower 
maxeig=double(real(eigenval(1))); %take the max real part of the first element of the array (which is the maxeig)

tol=10^-5;

if maxeig<-tol % max eig is negative 
    EIG=1; % Stable (S)
elseif abs(maxeig)<=tol % max eig is zero
    if ( abs(double(real(eigenval(2))))<=tol )
        if ( (abs(eigenval(2)-eigenval(3))<=tol) || (abs(eigenval(2)-eigenval(4))<=tol) || (abs(eigenval(2)-eigenval(5))<=tol) )% check if 2nd and 3rd/4th/5th complex eigv are equal - alg mult 2
            EIG=0;
            [V,D]=eigs(A2,5,'largestreal');
            if ( max( max(abs(real(V(:,2)-V(:,3)))), max(abs(imag(V(:,2)-V(:,3)))) )<=tol || ...
                 max( max(abs(real(V(:,2)-V(:,4)))), max(abs(imag(V(:,2)-V(:,4)))) )<=tol || ...
                 max( max(abs(real(V(:,2)-V(:,5)))), max(abs(imag(V(:,2)-V(:,5)))) )<=tol )
                disp('Unstable! --> geom. mult < alg mult !! ');
                EIG = -1; % Unstable since alg mult = 2 for complex eigval
            end
        else
            EIG=0; % MS if alg mult = 1 
        end
    else
        EIG=0; % Marginally stable (MS)
    end
else % max eig is positive 
    EIG=-1; % Unstable (U)
end

end