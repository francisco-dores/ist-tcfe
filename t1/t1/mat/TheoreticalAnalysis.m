close all
clear all


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


%%-----------------> MESH METHOD <-------------------

%%Coeficients matrix
A = [R1+R3+R4, -R4, -R3, 0; -R4, R4+R6+R7-Kc, 0, 0; R3*Kb, 0, 1-Kb*R3, 0; 0, 0, 0, 1];

%%Solution Matrix
b = [-Va; -V0; 0; Id];

%%Computing currents vector
I = inv(A)*b;


%%Show Solution

output_precision(10)

I1 = I(1);
I2 = I(2);
I3 = I(3);
I4 = I(4);

Ib=I3;
IR1=-I1;
IR2=I3;
IR3=I3-I1;
IR4=I2-I1;
IR5=I3-I4;
IR6=I2;
IR7=I2;



%%-----------------> NODE METHOD <-------------------

C=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,0,0,0,G3; ...
   0,G2+Kb,-G2,0,0,0,-Kb; ...
   0,-Kb,0,-G5,0,0,G5+Kb; ...
   0,0,0,0,-1,Kc*G6,1; ...
   0,0,0,0,G7,-G6-G7,0; ...
   0,G3,0,G5,0,-G6,-G3-G4-G5];
d=[Va;0;0;-Id;0;0;Id];

V=inv(C)*d;

printf('op_TAB\n');
printf('$I_b$ = %e\n$I_d$ = %e\n$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',Ib,Id,IR1,IR2,IR3,IR4,IR5,IR6,IR7);
printf('$V_1$ = %f\n$V_2$ = %f\n$V_3$ = %f\n$V_4$ = %f\n$V_5$ = %f\n$V_6$ = %f\n$V_7$ = %f\n', V(1),V(2),V(3),V(4),V(5),V(6),V(7));
printf('op_END');
