function E = equationmap_Tdip(y,par)
% This function is used to check if the equations are all satisfied
%global  r kappa  a ps kappap kappav  kappas
r = par(1); kappa = par(2); a = par(3); ps = par(4); 
kappap = par(5); kappav = par(6); kappas = par(7);

y(3)= y(3)+(1e-15)*(y(3)==0); % To avoid div by zero error we set lambda to very small num instead of 0

zita1=y(1); zita2=y(2); lambda=y(3); thetaf1=y(4); thetaf2=y(5);

%% dksi1/dt
ksi1dot=sin(thetaf1)+(1/24).*exp(1).^((sqrt(-1)*(-1)).*thetaf2).*pi.^2.*r.^2.*(( ...
  sqrt(-1)*(-3)).*ps.*(csch((1/2).*pi.*(lambda+sqrt(-1).*(zita1+(-1).*zita2))) ...
  .^2+sech((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+zita2))).^2)+(sqrt(-1)*3) ...
  .*exp(1).^((sqrt(-1)*2).*thetaf2).*ps.*(csch((1/2).*pi.*(lambda+(sqrt(-1)*( ...
  -1)).*(zita1+(-1).*zita2))).^2+sech((1/2).*pi.*(lambda+sqrt(-1).*(zita1+zita2))) ...
  .^2)+(-2).*exp(1).^(sqrt(-1).*thetaf2).*((-1)+3.*sec(pi.*zita1).^2).*sin( ...
  thetaf1));
%% dksi2/dt
ksi2dot=(sqrt(-1)*(-1/8)).*exp(1).^((sqrt(-1)*(-1)).*thetaf1).*pi.^2.*r.^2.*( csch((1/2).*pi.*...
    (lambda+sqrt(-1).*(zita1+(-1).*zita2))).^2+(-1).*exp(1).^(( ...
  sqrt(-1)*2).*thetaf1).*(csch((1/2).*pi.*(lambda+(sqrt(-...
1)*(-1)).*zita1+sqrt( ...
  -1).*zita2)).^2+sech((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*...
(zita1+zita2))).^2)+ ...
  sech((1/2).*pi.*(lambda+sqrt(-1).*(zita1+zita2))).^2)+...
  (1/12).*ps.*(12+pi.^2.* ...
  r.^2+(-3).*pi.^2.*r.^2.*sec(pi.*zita2).^2).*sin(thetaf2);
%% lambdadot
lambdadot=(1/24).*(96.*a.*zita1.^2+(-96).*a.*zita2.^2+(-24).*cos(thetaf1)+2.*pi.^2.* ...
  r.^2.*cos(thetaf1)+24.*ps.*cos(thetaf2)+(-2).*pi.^2.*r.^2.*ps.*cos(thetaf2)+6.* ...
  pi.^2.*r.^2.*cos(thetaf1).*sec(pi.*zita1).^2+(-6).*pi.^2.*r.^2.*ps.*cos( ...
  thetaf2).*sec(pi.*zita2).^2+(-3).*pi.^2.*r.^2.*cos(thetaf1).*sech((1/2).*pi.* ...
  (lambda+(sqrt(-1)*(-1)).*(zita1+zita2))).^2+3.*pi.^2.*r.^2.*ps.*cos(thetaf2).* ...
  sech((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+zita2))).^2+(-3).*pi.^2.* ...
  r.^2.*cos(thetaf1).*sech((1/2).*pi.*(lambda+sqrt(-1).*(zita1+zita2))).^2+3.* ...
  pi.^2.*r.^2.*ps.*cos(thetaf2).*sech((1/2).*pi.*(lambda+sqrt(-1).*(zita1+zita2))) ...
  .^2+(sqrt(-1)*(-3)).*pi.^2.*r.^2.*sech((1/2).*pi.*(lambda+(sqrt(-1)*( ...
  -1)).*(zita1+zita2))).^2.*sin(thetaf1)+(sqrt(-1)*3).*pi.^2.*r.^2.*sech((1/2) ...
  .*pi.*(lambda+sqrt(-1).*(zita1+zita2))).^2.*sin(thetaf1)+(sqrt(-1)*(-3)).*pi.^2.* ...
  r.^2.*ps.*sech((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+zita2))).^2.*sin( ...
  thetaf2)+(sqrt(-1)*3).*pi.^2.*r.^2.*ps.*sech((1/2).*pi.*(lambda+sqrt(-1).*( ...
  zita1+zita2))).^2.*sin(thetaf2)+3.*pi.^2.*r.^2.*csch((1/2).*pi.*(lambda+sqrt(-1) ...
  .*(zita1+(-1).*zita2))).^2.*(cos(thetaf1)+(-1).*ps.*cos(thetaf2)+(sqrt(-1)*(-1)) ...
  .*(sin(thetaf1)+(-1).*ps.*sin(thetaf2)))+3.*pi.^2.*r.^2.*csch((1/2).*pi.*( ...
  lambda+(sqrt(-1)*(-1)).*(zita1+(-1).*zita2))).^2.*(cos(thetaf1)+(-1).*ps.*cos(thetaf2) ...
  +sqrt(-1).*(sin(thetaf1)+(-1).*ps.*sin(thetaf2))));
%% T-DIPOLE
%thetaf2dot
theta1dot=(1/8).*exp(1).^((sqrt(-1)*(-1)).*(thetaf1+thetaf2)).*(64.*exp(1).^(sqrt( ...
  -1).*(thetaf1+thetaf2)).*a.*kappa.*zita1+64.*exp(1).^(sqrt(-1).*(thetaf1+thetaf2)).*a.* ...
  zita1.*cos(thetaf1).^2+(sqrt(-1)*(-1)).*exp(1).^((sqrt(-1)*3).*thetaf1+(sqrt( ...
  -1)*2).*thetaf2).*pi.^3.*r.^2.*ps.*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).* ...
  (zita1+(-1).*zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+(-1).*zita2) ...
  )).^2+sqrt(-1).*exp(1).^((sqrt(-1)*(-1)).*thetaf1).*pi.^3.*r.^2.*ps.* ...
  coth((1/2).*pi.*(lambda+sqrt(-1).*(zita1+(-1).*zita2))).*csch((1/2).*pi.*(lambda+ ...
  sqrt(-1).*(zita1+(-1).*zita2))).^2+(sqrt(-1)*(-1)).*exp(1).^(sqrt(-1).* ...
  thetaf1).*pi.^3.*r.^2.*ps.*cos(thetaf1).^2.*coth((1/2).*pi.*(lambda+(sqrt(-1)*( ...
  -1)).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))) ...
  .^2+sqrt(-1).*exp(1).^(sqrt(-1).*(thetaf1+2.*thetaf2)).*pi.^3.*r.^2.*ps.* ...
  cos(thetaf1).^2.*coth((1/2).*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2) ...
  .*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).^2+sqrt(-1).*exp(1).^(sqrt(-1).* ...
  thetaf1).*pi.^3.*r.^2.*ps.*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+ ...
  zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))).^2.*sin( ...
  thetaf1).^2+(sqrt(-1)*(-1)).*exp(1).^(sqrt(-1).*(thetaf1+2.*thetaf2)).*pi.^3.* ...
  r.^2.*ps.*coth((1/2).*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2).* ...
  pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).^2.*sin(thetaf1).^2+exp(1).^(sqrt(-1).* ...
  thetaf1).*pi.^3.*r.^2.*ps.*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+ ...
  zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))).^2.*sin(2.* ...
  thetaf1)+exp(1).^(sqrt(-1).*(thetaf1+2.*thetaf2)).*pi.^3.*r.^2.*ps.*coth((1/2) ...
  .*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+sqrt(-1).*(1+ ...
  zita1+zita2))).^2.*sin(2.*thetaf1)+(1/2).*exp(1).^(sqrt(-1).*thetaf2).*pi.^3.* ...
  r.^2.*sec(pi.*zita1).^2.*tan(pi.*zita1)+(-1/2).*exp(1).^(sqrt(-1).*(4.* ...
 thetaf1+thetaf2)).*pi.^3.*r.^2.*sec(pi.*zita1).^2.*tan(pi.*zita1)+exp(1).^(sqrt( ...
  -1).*thetaf2).*pi.^3.*r.^2.*cos(thetaf1).^2.*sec(pi.*zita1).^2.*tan(pi.*zita1)+ ...
  exp(1).^(sqrt(-1).*(2.*thetaf1+thetaf2)).*pi.^3.*r.^2.*cos(thetaf1).^2.*sec( ...
  pi.*zita1).^2.*tan(pi.*zita1)+(sqrt(-1)*2).*exp(1).^(sqrt(-1).*thetaf2).* ...
  pi.^3.*r.^2.*cos(thetaf1).*sec(pi.*zita1).^2.*sin(thetaf1).*tan(pi.*zita1)+(-1) ...
  .*exp(1).^(sqrt(-1).*thetaf2).*pi.^3.*r.^2.*sec(pi.*zita1).^2.*sin(thetaf1) ...
  .^2.*tan(pi.*zita1)+(-1).*exp(1).^(sqrt(-1).*(2.*thetaf1+thetaf2)).*pi.^3.* ...
  r.^2.*sec(pi.*zita1).^2.*sin(thetaf1).^2.*tan(pi.*zita1))+kappas*sin(thetaf1)+(1+(-1).*cos(thetaf2+atan(lambda.^(-1).*(zita1+(-1).*zita2)))).*(kappav.*ps.*sin(thetaf1+( ...
  -1).*thetaf2)+kappap.*(lambda.^2+(zita1+(-1).*zita2).^2).^(1/2).*sin(thetaf2+atan(lambda.^(-1) ...
  .*(zita1+(-1).*zita2))));

% thetaf2dot
theta2dot=(1/8).*exp(1).^((sqrt(-1)*(-1)).*(thetaf1+thetaf2)).*(64.*exp(1).^(sqrt( ...
  -1).*(thetaf1+thetaf2)).*a.*kappa.*zita2+64.*exp(1).^(sqrt(-1).*(thetaf1+thetaf2)).*a.* ...
  zita2.*cos(thetaf2).^2+sqrt(-1).*exp(1).^((sqrt(-1)*2).*thetaf1+(sqrt(-1)*3) ...
  .*thetaf2).*pi.^3.*r.^2.*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+(-1) ...
  .*zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(zita1+(-1).*zita2))).^2+( ...
  sqrt(-1)*(-1)).*exp(1).^((sqrt(-1)*(-1)).*thetaf2).*pi.^3.*r.^2.*coth( ...
  (1/2).*pi.*(lambda+sqrt(-1).*(zita1+(-1).*zita2))).*csch((1/2).*pi.*(lambda+sqrt( ...
  -1).*(zita1+(-1).*zita2))).^2+(sqrt(-1)*(-1)).*exp(1).^(sqrt(-1).*(2.* ...
  thetaf1+thetaf2)).*pi.^3.*r.^2.*cos(thetaf2).^2.*coth((1/2).*pi.*(lambda+(sqrt(-1)* ...
  (-1)).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2)) ...
  ).^2+sqrt(-1).*exp(1).^(sqrt(-1).*thetaf2).*pi.^3.*r.^2.*cos(thetaf2).^2.* ...
  coth((1/2).*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+sqrt( ...
  -1).*(1+zita1+zita2))).^2+(-2).*exp(1).^(sqrt(-1).*(2.*thetaf1+thetaf2)).* ...
  pi.^3.*r.^2.*cos(thetaf2).*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+ ...
  zita2))).*csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))).^2.*sin( ...
  thetaf2)+(-2).*exp(1).^(sqrt(-1).*thetaf2).*pi.^3.*r.^2.*cos(thetaf2).*coth(( ...
  1/2).*pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+sqrt(-1).*( ...
  1+zita1+zita2))).^2.*sin(thetaf2)+sqrt(-1).*exp(1).^(sqrt(-1).*(2.*thetaf1+thetaf2)) ...
  .*pi.^3.*r.^2.*coth((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))).* ...
  csch((1/2).*pi.*(lambda+(sqrt(-1)*(-1)).*(1+zita1+zita2))).^2.*sin(thetaf2).^2+( ...
  sqrt(-1)*(-1)).*exp(1).^(sqrt(-1).*thetaf2).*pi.^3.*r.^2.*coth((1/2).* ...
  pi.*(lambda+sqrt(-1).*(1+zita1+zita2))).*csch((1/2).*pi.*(lambda+sqrt(-1).*(1+zita1+ ...
  zita2))).^2.*sin(thetaf2).^2+(1/2).*exp(1).^(sqrt(-1).*thetaf1).*pi.^3.* ...
  r.^2.*ps.*sec(pi.*zita2).^2.*tan(pi.*zita2)+(-1/2).*exp(1).^(sqrt(-1).*( ...
  thetaf1+4.*thetaf2)).*pi.^3.*r.^2.*ps.*sec(pi.*zita2).^2.*tan(pi.*zita2)+exp(1) ...
  .^(sqrt(-1).*thetaf1).*pi.^3.*r.^2.*ps.*cos(thetaf2).^2.*sec(pi.*zita2).^2.* ...
  tan(pi.*zita2)+exp(1).^(sqrt(-1).*(thetaf1+2.*thetaf2)).*pi.^3.*r.^2.*ps.*cos( ...
  thetaf2).^2.*sec(pi.*zita2).^2.*tan(pi.*zita2)+(sqrt(-1)*2).*exp(1).^(sqrt( ...
  -1).*thetaf1).*pi.^3.*r.^2.*ps.*cos(thetaf2).*sec(pi.*zita2).^2.*sin(thetaf2).* ...
  tan(pi.*zita2)+(-1).*exp(1).^(sqrt(-1).*thetaf1).*pi.^3.*r.^2.*ps.*sec( ...
  pi.*zita2).^2.*sin(thetaf2).^2.*tan(pi.*zita2)+(-1).*exp(1).^(sqrt(-1).*( ...
  thetaf1+2.*thetaf2)).*pi.^3.*r.^2.*ps.*sec(pi.*zita2).^2.*sin(thetaf2).^2.*tan( ...
  pi.*zita2))+kappas*ps*sin(thetaf2)+(1+(-1).*cos(thetaf2+atan(lambda.^(-1).*(zita1+(-1).*zita2)))).*(kappav.*ps.*sin(thetaf1+( ...
  -1).*thetaf2)+kappap.*(lambda.^2+(zita1+(-1).*zita2).^2).^(1/2).*sin(thetaf2+atan(lambda.^(-1) ...
  .*(zita1+(-1).*zita2))));

tol=10^-5; %tolerance can be set

if max([abs(ksi1dot) abs(ksi2dot) abs(lambdadot) abs(theta1dot) abs(theta2dot)])<=tol
% if abs(ksi1eq)<=tol && abs(ksi2eq)<=tol && abs(lambdaeq)<=tol && abs(theta1eq)<=tol && abs(theta2eq)<=tol
    E=1;
else 
    E=0;
end

end