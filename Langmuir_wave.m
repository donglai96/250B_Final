clear;clc;
zeta = @(x) faddeeva(x)*1i*sqrt(pi);
f = @(x,k)1 + (k*k) + x*zeta(x);w=[];
kmin = 0.1; dk1 = 0.1; kmid = 1; dk2 = 1; kmax = 10.0;
k = [kmin:dk1:kmid,(kmid + dk2):dk2:kmax];
c = 1-0.1i;
for kk = k
    options = optimset('Display','off');
    %x = fsolve(f, 1-0.1i,options,kk)*sqrt(2)*kk;
    
    c = fsolve(f, c, options, kk);
    x = c*sqrt(2)*kk;
    w=[w,x];
end
wre = real(w);wie = imag(w);
wrt = 1.0 + 1.5.*k.*k;
wit = -sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3);

loglog(k,wre,'r-',k,-wie,'r--');
legend('\omega-r','-\gamma','Location','SouthEast');
xlabel('k\lambda_D');ylabel('\omega/\omega_p');
title('Langmuir Dispersion Relation')
xlim([kmin,kmax]);
ylim([0.0001,100]);
