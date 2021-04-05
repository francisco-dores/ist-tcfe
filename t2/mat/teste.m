close all
clear all

format long

%%Variable Values
R1 = 1.02597459645e3; 
R2 = 2.02178008702e3; 
R3 = 3.03144887426e3; 
R4 = 4.13109833342e3; 
R5 = 3.09601431108e3; 
R6 = 2.01920057699e3; 
R7 = 1.02918842978e3; 
Vs = 5.1256272592; 
C = 1.011814928e-6; 
Kb = 7.28538907285e-3; 
Kd = 8.0919603219e3;


G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;


%%-----------------> NODE METHOD to calculate voltages t<0 <-------------------

C=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,G3,0,0,0; ...
   0,G2+Kb,-G2,-Kb,0,0,0; ...
   0,G3,0,-G3-G4-G5,G5,G7,-G7; ...
   0,-Kb,0,Kb+G5,-G5,0,0; ...
   0,0,0,0,0,-G6-G7,G7; ...
   0,0,0,1,0,Kd*G6,-1];
d=[Vs;0;0;0;0;0;0];

V=inv(C)*d;

V1=V(1);
V2=V(2);
V3=V(3);
V4=0;
V5=V(4);
V6=V(5);
V7=V(6);
V8=V(7);


%printf('$I_b$ = %e\n$I_d$ = %e\n$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',Ib,Id,IR1,IR2,IR3,IR4,IR5,IR6,IR7);


printf('op_TAB\n');
printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6,V7,V8);


%%-----------------> Calculate equivalent resistor <-------------------

Vx=V6-V8;

%C=[1,0,0,0,0,0,0; ...
%   G1,-G1-G2-G3,G2,G3,0,0,0; ...
%   0,G2+Kb,-G2,-Kb,0,0,0; ...
%   0,0,0,1,0,Kd*G6,-1; ...
%   0,0,0,0,1,0,-1; ...
%   0,0,0,0,0,-G6-G7,G7; ...
%   G1,-G1,0,-G4,0,-G6,0];
Vs=0;  
C=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,G3,0,0,0; ...
   0,G2+Kb,-G2,-Kb,0,0,0; ...
   0,0,0,1,0,Kd*G6,-1; ...
   0,0,0,0,1,0,-1; ...
   0,0,0,0,0,-G6-G7,G7; ...
   G4,G3,0,-G3-G4,0,G6+G7,-G7];   

   
d=[Vs;0;0;0;Vx;0;0];

V=inv(C)*d;

V1=V(1);
V2=V(2);
V3=V(3);
V4=0;
V5=V(4);
V6=V(5);
V7=V(6);
V8=V(7);


Ix=G5*(V5-V6)-Kb*(V2-V5);
%Ix=-G7*(V7-V8)+Kd*G6*V7;
Req=Vx/Ix;


%printf('$I_b$ = %e\n$I_d$ = %e\n$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',Ib,Id,IR1,IR2,IR3,IR4,IR5,IR6,IR7);


printf('op_TAB\n');
printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6,V7,V8);
printf('op_END');

%%-----------------> Calculate natural solution <-------------------

t=0:1e-6:20e-3;

wn=-1/Req*C

vc_n=Vx*exp(-wn*t);

nat_sol = figure ();
plot (t*1000, vc_n, "g"); 
xlabel ("t [ms]");
ylabel ("V_(6n) [V]");
print (nat_sol, "nat_sol.eps", "-depsc");


%%-----------------> Calculate forced solution <-------------------

%f=1000;
%w=2*pi*f; %certo?
%Zc=1/(jwC);

%CGain=Zc/(Zc+Req)*exp(-pi/2*j);
%Gain=abs(CGain);
%Phase=angle(CGain);

sin(w*t)=exp(jwt)*exp(-pi/2*j)




%vc_f=Gain*cos(w*t+Phase);






