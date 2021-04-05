close all
clear all

%%EXAMPLE SYMBOLIC COMPUTATIONS

pkg load symbolic

syms t
syms R
syms C
syms vs(t)
syms vc(t)
syms vc_n(t) %natural solution
syms i(t)
syms A
syms wn


R=1e3 %Ohm
C=100e-9 %F
A = V6-V8

i(t)=C*diff(vc,t)

printf("\n\nKVL equation:\n");

vs(t) = R*i(t)+vc(t)

printf("\n\nSolution is of the form:\n");

vc(t) = vc_n(t) + vc_f(t)

printf("\n\nNatural solution:\n");

vs(t) = 0 %no excitation
i_n(t) = C*diff(vc_n, t)


printf("\n\n Natural solution is of the form:\n");
vc_n(t) = A*exp(wn*t)

R*i_n(t)+vc_n(t) == 0

R*C*wn*vc_n(t)+vc_n(t) == 0

R*C*wn+1==0

wn=-1/R*C


%time axis: 0 to 20ms with 1us steps
t=0:1e-6:20e-3; %s

hf = figure ();
plot (t*1000, vs, "g");
hold on;
plot (t*1000, vcn, "b");

xlabel ("t[ms]");
ylabel ("vs(t), vcn(t) [V]");
print (hf, "forced.eps", "-depsc");
