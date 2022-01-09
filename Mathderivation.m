%% Symbolic formulation of lognormal distribution
clear
clc

%% Lognormal as symbolic function
syms f a b x
f(x) = (1/(x*b*sqrt(2*sym(pi))))*exp((-(log(x)-a)^2)/(2*b^2));
pretty(f)

%% Lognormal as anonymous function
aa = 3;
bb = 0.3;
ff = @(xx) (sqrt(2).*exp(-(aa - log(xx)).^2./(2.*bb.^2)))./(2.*bb.*xx*sqrt(pi));

%% Derivative with respect to a
dfa = diff(f,a);
pretty(dfa)

%% Derivative with respect to b
dfb = diff(f,b);
pretty(dfb)

%% Prepare input file
clf
aa = 2.5;
bb = 0.254;
ff = @(xx) (sqrt(2).*exp(-(aa - log(xx)).^2./(2.*bb.^2)))./(2.*bb.*xx*sqrt(pi));
plot(ff(1:100))
%%
fid = fopen('urf01.dat','w');
fprintf(fid,'%.10f\n', ff(1:100));
fclose(fid);