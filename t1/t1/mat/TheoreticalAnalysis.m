close all
clear all

%%-----------------> MESH METHOD <-------------------
%%Variable Values
R1 = 1.02597459645e3;
R2 = 2.02178008702e3;
R3 = 3.03144887426e3;
R4 = 4.13109833342e3;
R5 = 3.09601431108e3;
R6 = 2.01920057699e3;
R7 = 1.02918842978e3;
Va = 5.1256272592;
V0 = 0;
Kb = 7.28538907285e-3;
Kc = 8.0919603219e3;
Id = 1.011814928e-3;

G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

printf("\n\n")


%%Coeficients matrix
A = [R1+R3+R4, -R4, -R3, 0; -R4, R4+R6+R7-Kc, 0, 0; R3*Kb, 0, 1-Kb*R3, 0; 0, 0, 0, 1];

%%Solution Matrix
b = [-Va; -V0; 0; Id];

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

V2= Va+R1*I1
V3= V2+R2*I3

Vb = -I1*R3
Ib = Kb*Vb

printf("\n")

Ic = I2
Vc = Kc*Ic

printf("\n\n")

%%-----------------> NODE METHOD <-------------------

C=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,0,0,0,G3; ...
   0,G2+Kb*G3,-G2,0,0,0,-Kb*G3; ...
   0,-Kb*G3,0,-G5,0,0,G5+Kb*G3; ...
   0,0,0,0,1,Kc*G6,-1; ...
   0,0,0,0,G7,-G6-G7,0; ...
   0,G3,0,G5,-G7,G7,-G3-G4];
d=[Va;0;0;-Id;0;0;Id];

V=inv(C)*d

printf("\n\n")

Vb = V(2)-V(7)

printf("\n")

Vc = V(7)-V(5)

%%Plot

hf = figure ();
plot (t*1000, vi, "g");
hold on;
plot (t*1000, vo, "b");

xlabel ("t[ms]");
ylabel ("vi(t), vo(t) [V]");
print (hf, "forced.eps", "-depsc");










