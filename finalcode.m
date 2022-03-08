clc
clear all 
close all 

%% Parameters of Blasius Equation

Pr = input('Enter the value of prandtl number(Enter between 10e-3 to 10e3 for quick solution): ');


%% Shooting method to find out f11 at eta equal to zero
fd = @(x, f, f1, f11) f1;
fdd = @(x, f, f1, f11) f11;
ftd = @(x, f, f1, f11) -f*f11;
h=0.05;
eta = 0:h:100;
x = 0:h:100;
f(1) = 0;
f1(1) = 0;
j=1;
f11g(1)=1.5;
f11gp = 2;
f11gn = 0.1;
while f11gp-f11gn>0.00001
    f11(1)=f11g(j);
for i = 1:(length(eta)-1)
a = h.*[fd(eta(i), f(i), f1(i), f11(i)), fdd(eta(i), f(i), f1(i), f11(i)), ftd(eta(i), f(i), f1(i), f11(i))];
b = h.*[fd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), fdd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), ftd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2)];
c = h.*[fd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), fdd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), ftd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2)];
d = h.*[fd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), fdd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), ftd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3))];
f11(i+1) = f11(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
f1(i+1) = f1(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
f(i+1) = f(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end
if (f1(1 + 10/h) < 1)
    f11gn= f11g(j);
elseif (f1(1 + 10/h) > 1)
    f11gp = f11g(j);
end
    f11g(j+1) = (f11gn+f11gp)/2;
    j=j+1;
end 
disp(f11g(j))
f11(1)=f11g(j);
%% Shooting method to get phi1
etamax= 10/sqrt(Pr);
h = etamax/10000;
eta = 0:h:etamax;
x = 0:h:etamax;
phid= @(x,f,f1,phi,phi1) phi1;
phidd= @(x,f,f1,phi,phi1) -Pr*(f*phi1 - phi*f1);
f(1) = 0;
f1(1) = 0;
phi1(1) = -1;
phi1g(1)= 28;
phi1gp = 90;
phi1gn = 0.1;
j=1;
while phi1gp-phi1gn>0.00001
    phi(1)=phi1g(j);
for i = 1:(length(eta)-1)
a = h.*[fd(eta(i), f(i), f1(i), f11(i)), fdd(eta(i), f(i), f1(i), f11(i)), ftd(eta(i), f(i), f1(i), f11(i))];
b = h.*[fd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), fdd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), ftd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2)];
c = h.*[fd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), fdd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), ftd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2)];
d = h.*[fd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), fdd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), ftd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3))];
e = h.*[phid(eta(i), f(i), f1(i), phi(i), phi1(i)), phidd(eta(i), f(i), f1(i), phi(i), phi1(i))];
ff= h.*[phid(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, phi(i)+e(1)/2, phi1(i)+e(2)/2), phidd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, phi(i)+e(1)/2, phi1(i)+e(2)/2)];
g=  h.*[phid(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, phi(i)+ff(1)/2, phi1(i)+ff(2)/2), phidd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, phi(i)+ff(1)/2, phi1(i)+ff(2)/2)];
hh= h.*[phid(eta(i)+h, f(i)+c(1), f1(i)+c(2), phi(i)+g(1), phi1(i)+g(2)), phidd(eta(i)+h, f(i)+c(1), f1(i)+c(2), phi(i)+g(1), phi1(i)+g(2))];
f11(i+1) = f11(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
f1(i+1) = f1(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
f(i+1) = f(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
phi1(i+1) = phi1(i) + 1/6*(e(2)+2*ff(2)+2*g(2)+hh(2));
phi(i+1) = phi(i) + 1/6*(e(1)+2*ff(1)+2*g(1)+hh(1));
end
if (phi(1 + etamax/h) < 0)
    phi1gn= phi1g(j);
elseif (phi(1 + etamax/h) > 0)
    phi1gp = phi1g(j);
end
    phi1g(j+1) = (phi1gn+phi1gp)/2;
    j=j+1;
end 
disp(phi1g(j));
phi(1)=phi1g(j);
%% Numerical Solution of Blasius Equation Using Runge-Kutta
%etamax= 0.2/Pr;
etamax= 10/sqrt(Pr);
h = etamax/10000;
eta = 0:h:etamax;
x = 0:h:etamax;
fd = @(x, f, f1, f11) f1;
fdd = @(x, f, f1, f11) f11;
ftd = @(x, f, f1, f11) -f*f11;
phid= @(x,f,f1,phi,phi1) phi1;
phidd= @(x,f,f1,phi,phi1) -Pr*(f*phi1 - phi*f1);
f(1) = 0;
f1(1) = 0;
phi1(1) = -1;
for i = 1:(length(eta)-1)
a = h.*[fd(eta(i), f(i), f1(i), f11(i)), fdd(eta(i), f(i), f1(i), f11(i)), ftd(eta(i), f(i), f1(i), f11(i))];
b = h.*[fd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), fdd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2), ftd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, f11(i)+a(3)/2)];
c = h.*[fd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), fdd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2), ftd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, f11(i)+b(3)/2)];
d = h.*[fd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), fdd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3)), ftd(eta(i)+h, f(i)+c(1), f1(i)+c(2), f11(i)+c(3))];
e = h.*[phid(eta(i), f(i), f1(i), phi(i), phi1(i)), phidd(eta(i), f(i), f1(i), phi(i), phi1(i))];
ff= h.*[phid(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, phi(i)+e(1)/2, phi1(i)+e(2)/2), phidd(eta(i)+h/2, f(i)+a(1)/2, f1(i)+a(2)/2, phi(i)+e(1)/2, phi1(i)+e(2)/2)];
g=  h.*[phid(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, phi(i)+ff(1)/2, phi1(i)+ff(2)/2), phidd(eta(i)+h/2, f(i)+b(1)/2, f1(i)+b(2)/2, phi(i)+ff(1)/2, phi1(i)+ff(2)/2)];
hh= h.*[phid(eta(i)+h, f(i)+c(1), f1(i)+c(2), phi(i)+g(1), phi1(i)+g(2)), phidd(eta(i)+h, f(i)+c(1), f1(i)+c(2), phi(i)+g(1), phi1(i)+g(2))];
f11(i+1) = f11(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
f1(i+1) = f1(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
f(i+1) = f(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
phi1(i+1) = phi1(i) + 1/6*(e(2)+2*ff(2)+2*g(2)+hh(2));
phi(i+1) = phi(i) + 1/6*(e(1)+2*ff(1)+2*g(1)+hh(1));
end

%% Plotting and Visualization
plot( eta, f1, 'LineWidth', 2)
xlim([0 etamax])
xlabel('eta')
ylabel('phi')
title('Plot of phi for low Pr')
grid on


