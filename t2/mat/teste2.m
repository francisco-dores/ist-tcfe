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

H=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,G3,0,0,0; ...
   0,G2+Kb,-G2,-Kb,0,0,0; ...
   0,G3,0,-G3-G4-G5,G5,G7,-G7; ...
   0,-Kb,0,Kb+G5,-G5,0,0; ...
   0,0,0,0,0,-G6-G7,G7; ...
   0,0,0,1,0,Kd*G6,-1];
d=[Vs;0;0;0;0;0;0];

V=inv(H)*d;

V1=V(1);
V2=V(2);
V3=V(3);
V4=0;
V5=V(4);
V6_b=V(5);
V7=V(6);
V8=V(7);


%printf('$I_b$ = %e\n$I_d$ = %e\n$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',Ib,Id,IR1,IR2,IR3,IR4,IR5,IR6,IR7);


printf('op_TAB\n');
printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6_b,V7,V8);
printf('op_END\n');

%%-----------------> Calculate equivalent resistor <-------------------

Vx=V6_b-V8;

%H=[1,0,0,0,0,0,0; ...
%   G1,-G1-G2-G3,G2,G3,0,0,0; ...
%   0,G2+Kb,-G2,-Kb,0,0,0; ...
%   0,0,0,1,0,Kd*G6,-1; ...
%   0,0,0,0,1,0,-1; ...
%   0,0,0,0,0,-G6-G7,G7; ...
%   G1,-G1,0,-G4,0,-G6,0];
Vs=0;  
H=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,G3,0,0,0; ...
   0,G2+Kb,-G2,-Kb,0,0,0; ...
   0,0,0,1,0,Kd*G6,-1; ...
   0,0,0,0,1,0,-1; ...
   0,0,0,0,0,-G6-G7,G7; ...
   G4,G3,0,-G3-G4,0,G6+G7,-G7];   

   
d=[Vs;0;0;0;Vx;0;0];

V=inv(H)*d;

V1=V(1);
V2=V(2);
V3=V(3);
V4=0;
V5=V(4);
V6_0=V(5);
V7=V(6);
V8=V(7);


Ix=G5*(V5-V6_0)-Kb*(V2-V5);
%Ix=-G7*(V7-V8)+Kd*G6*V7;
Req=Vx/Ix; 
Req=-Req; %Não sabemos ainda o erro do Ix


%printf('$I_b$ = %e\n$I_d$ = %e\n$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',Ib,Id,IR1,IR2,IR3,IR4,IR5,IR6,IR7);


printf('op_TAB\n');
printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6_0,V7,V8);
printf('op_END\n');

%%-----------------> Calculate natural solution <-------------------

t=0:1e-6:20e-3;

Vx=V6_0-V8;

wn=-1/(Req*C);

v6_n=Vx*exp(wn*t);

nat_sol = figure ();
plot (t*1000, v6_n, "g"); 
xlabel ("t [ms]");
ylabel ("V_(6n) [V]");
print (nat_sol, "nat_sol.eps", "-depsc");


%%-----------------> Calculate forced solution <-------------------
f=1000;
w=2*pi*f;
Zc=1/(j*w*C);
Gc=1/Zc;
Vs=exp(-pi/2*j);


H=[1,0,0,0,0,0,0; ...
   G1,-G1-G2-G3,G2,G3,0,0,0; ...
   0,G2+Kb,-G2,-Kb,0,0,0; ...
   0,G3,0,-G3-G4-G5,G5+Gc,G7,-Gc-G7; ...
   0,-Kb,0,Kb+G5,-G5-Gc,0,Gc; ...
   0,0,0,0,0,-G6-G7,G7; ...
   0,0,0,1,0,Kd*G6,-1];
d=[Vs;0;0;0;0;0;0];

V=inv(H)*d;

V1=V(1);
V2=V(2);
V3=V(3);
V4=0;
V5=V(4);
V6_a=V(5);
V7=V(6);
V8=V(7);


printf('op_TAB\n');
printf('$V_1$ = %fe^{%fj}\n$V_2$ = %fe^{%fj}\n$V_3$ = %fe^{%fj}\n$V_4$ = %fe^{%fj}\n$V_5$ = %fe^{%fj}\n$V_6$ = %fe^{%fj}\n$V_7$ = %fe^{%fj}\n$V_8$ = %fe^{%fj}\n', abs(V1),angle(V1),abs(V2),angle(V2),abs(V3),angle(V3),abs(V4),angle(V4),abs(V5),angle(V5),...
abs(V6_a),angle(V6_a),abs(V7),angle(V7),abs(V8),angle(V8));
printf('op_END\n');

%%-----------------> Calculate total solution <-------------------


Gain=abs(V6_a);
Phase=angle(V6_a);

v6_forc=Gain*cos(w*t+Phase);

v6=v6_n+v6_forc;

vs=sin(w*t);

tb=-5e-3:1e-6:0;

V6_b=ones(1,length(tb))*V6_b;
Vs=ones(1,length(tb))*Vs;


tot_sol = figure ();
plot (t*1000, v6, "g");  
hold on;
plot (t*1000, vs, "b");
hold on;
plot (tb*1000, V6_b, "g"); 
%hold on;
%plot (0, V6_0, "g"); 
hold on;
plot (tb*1000, Vs, "b");
xlabel ("t [ms]");
ylabel ("Voltage [V]");
legend("v6","vs");
print (tot_sol, "tot_sol.eps", "-depsc");



%%-----------------> Frequency Response <-------------------

w=0.1*2*pi:100:1e6*2*pi;

%v6=(w.*C).*sqrt(1./(Req.^2.+w.^2.*C.^2));
%v6=1./(sqrt(1./(1+w.^2.*Req.^2.*C^2)))

freq_db = figure ();
%plot (log10(w/2/pi), 20*log10(abs(Vs)), "g");   %confirmar Vc????
%hold on;
%plot (log10(w/2/pi), 20*log10(abs(Vc)), "g");
%hold on;
plot (log10(w/2/pi), 20*log10(abs(v6)), "g");
xlabel ("log_10(f) [Hz]");
ylabel ("Magnitude v_s(f), v_c(f), v_6(f) [dB]");
print (freq_db, "freq_db.eps", "-depsc");


freq_ph = figure ();
plot (log10(w/2/pi), (ph(Vs)+pi/2)*180/pi, "g");   %confirmar Vc???? ph???
hold on;
plot (log10(w/2/pi), (ph(Vc)+pi/2)*180/pi, "g");
hold on;
plot (log10(w/2/pi), (ph(V6)+pi/2)*180/pi, "g");
xlabel ("log_10(f) [Hz]");
ylabel ("Phase v_s(f), v_c(f), v_6(f) [degrees]");
print (freq_ph, "freq_ph.eps", "-depsc");
















