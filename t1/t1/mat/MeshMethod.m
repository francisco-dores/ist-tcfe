close all
clear all

%%-----------------> MESH METHOD <-------------------

%%Variable Identification

pkg load symbolic

syms R1
syms R2
syms R3
syms R4
syms R5
syms R6
syms R7
syms Va
syms Vb
syms V0
syms Kb
syms Kc
syms Ic
syms Id
syms I1
syms I2
syms I3
syms I4


%%Variable Values
R1 = 1.02597459645;
R2 = 2.02178008702;
R3 = 3.03144887426;
R4 = 4.13109833342;
R5 = 3.09601431108;
R6 = 2.01920057699;
R7 = 1.02918842978;
Va = 5.1256272592;
V0 = 0;
Kb = 7.28538907285;
Kc = 8.0919603219;
Id = 1.011814928;

printf("\n\n")

%%Coeficients matrix
A = [R1+R3+R4, -R4, -R3, 0; -R4, R4+R6+R7-Kc, 0, 0; R1+R3, R6+R7-Kc, 0, -R3; 0, 0, 0, 1];

%%Solution Matrix
b = [-Va; -V0; -Va-V0; Id];

%%Inverting Matrix A
AI = inv(A);

printf("\n")

%%Computing currents vector
I = AI*b;


%%Show Solution

output_precision(10)

I1 = I(1)
I2 = I(2)
I3 = I(3)
I4 = I(4)

printf("\n\n")

Vb = -I1*R3
Ib = Kb*Vb

printf("\n")

Ic = I2
Vc = Kc*Ic

printf("\n\n")


%%Plot

hf = figure ();
plot (t*1000, vi, "g");
hold on;
plot (t*1000, vo, "b");

xlabel ("t[ms]");
ylabel ("vi(t), vo(t) [V]");
print (hf, "forced.eps", "-depsc");
