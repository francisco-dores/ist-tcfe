close all
clear all

pkg load symbolic

format long

values=dlmread("../data.txt");

%%Variable Values
R1 = values(3,4)*1000;
R2 = values(4,3)*1000;
R3 = values(5,3)*1000;
R4 = values(6,3)*1000;
R5 = values(7,3)*1000;
R6 = values(8,3)*1000;
R7 = values(9,3)*1000;
Vs = values(10,3);
C = values(11,3)*0.000001;
Kb = values(12,3)*0.001;
Kd = values(13,3)*1000; 

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

%printf('op_TAB\n');
%printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6_b,V7,V8);
%printf('op_END\n');

I1=G1*(V1-V2);
I2=G2*(V3-V2);
I3=G3*(V2-V5);
I4=G4*V5;
I5=G5*(V5-V6_b);
I6=-G6*V7;
I7=I6;
Ib=I2;
Ic=I5-Ib;
I_vd=I6+Ic;
I_vs=I4-I6;

tab=fopen("volt_tb0.tex", "w");
fprintf(tab, "$V_1$ & %.11f \\\\ \\hline \n", V1);
fprintf(tab, "$V_2$ & %.11f \\\\ \\hline \n", V2);
fprintf(tab, "$V_3$ & %.11f \\\\ \\hline \n", V3);
fprintf(tab, "$V_4$ & %.11f \\\\ \\hline \n", V4);
fprintf(tab, "$V_5$ & %.11f \\\\ \\hline \n", V5);
fprintf(tab, "$V_6$ & %.11f \\\\ \\hline \n", V6_b);
fprintf(tab, "$V_7$ & %.11f \\\\ \\hline \n", V7);
fprintf(tab, "$V_8$ & %.11f \\\\ \\hline \n", V8);
fclose(tab);

tab=fopen("curr_tb0.tex", "w");
fprintf(tab, "$I_1$ & %.11f \\\\ \\hline \n", I1);
fprintf(tab, "$I_2$ & %.11f \\\\ \\hline \n", I2);
fprintf(tab, "$I_3$ & %.11f \\\\ \\hline \n", I3);
fprintf(tab, "$I_4$ & %.11f \\\\ \\hline \n", I4);
fprintf(tab, "$I_5$ & %.11f \\\\ \\hline \n", I5);
fprintf(tab, "$I_6$ & %.11f \\\\ \\hline \n", I6);
fprintf(tab, "$I_7$ & %.11f \\\\ \\hline \n", I7);
fprintf(tab, "$I_b$ & %.11f \\\\ \\hline \n", Ib);
fprintf(tab, "$I_c$ & %.11f \\\\ \\hline \n", Ic);
fprintf(tab, "$I_{V_d}$ & %.11f \\\\ \\hline \n", I_vd);
fprintf(tab, "$I_{V_s}$ & %.11f \\\\ \\hline \n", I_vs);
fclose(tab);


%%-----------------> Calculate equivalent resistor <-------------------

Vx=V6_b-V8;

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

V61=V(5);
V81=V(7);

Ix=G5*(V5-V6_0)-Kb*(V2-V5);
%Ix=-G7*(V7-V8)+Kd*G6*V7;
Req=Vx/Ix; 
Req=-Req; %Não sabemos ainda o erro do Ix

printf('op_TAB\n');
%printf('$I_{R1}$ = %e\n$I_{R2}$ = %e\n$I_{R3}$ = %e\n$I_{R4}$ = %e\n$I_{R5}$ = %e\n$I_{R6}$ = %e\n$I_{R7}$ = %e\n',IR1,IR2,IR3,IR4,IR5,IR6,IR7);
printf('$V_1$ = %.11f\n$V_2$ = %.11f\n$V_3$ = %.11f\n$V_4$ = %.11f\n$V_5$ = %.11f\n$V_6$ = %.11f\n$V_7$ = %.11f\n$V_8$ = %.11f\n', V1,V2,V3,V4,V5,V6_0,V7,V8);
printf('op_END\n');

tab=fopen("volt_vs0.tex", "w");
fprintf(tab, "$V_1$ & %.11f \\\\ \\hline \n", V1);
fprintf(tab, "$V_2$ & %.11f \\\\ \\hline \n", V2);
fprintf(tab, "$V_3$ & %.11f \\\\ \\hline \n", V3);
fprintf(tab, "$V_4$ & %.11f \\\\ \\hline \n", V4);
fprintf(tab, "$V_5$ & %.11f \\\\ \\hline \n", V5);
fprintf(tab, "$V_6$ & %.11f \\\\ \\hline \n", V6_0);
fprintf(tab, "$V_7$ & %.11f \\\\ \\hline \n", V7);
fprintf(tab, "$V_8$ & %.11f \\\\ \\hline \n", V8);
fclose(tab);

tab=fopen("r_eq.tex", "w");
fprintf(tab, "$I_x$ & %.11f A \\\\ \\hline \n", Ix);
fprintf(tab, "$V_x$ & %.11f V \\\\ \\hline \n", Vx);
fprintf(tab, "$R_{eq}$ & %.11f Ohm \\\\ \\hline \n", Req);
fclose(tab);


%%-----------------> Calculate natural solution <-------------------

t=0:1e-6:20e-3;

Vx=V6_0-V8;

wn=-1/(Req*C);

v6_n=Vx*exp(wn*t);

nat_sol = figure ();
plot (t*1000, v6_n, "r"); 
xlabel ("t [ms]");
ylabel ("V_{6n} [V]");
title ("Natural Solution");
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


tab=fopen("v6f.tex", "w");
fprintf(tab, "$V_1$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V1),angle(V1));
fprintf(tab, "$V_2$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V2),angle(V2));
fprintf(tab, "$V_3$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V3),angle(V3));
fprintf(tab, "$V_4$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V4),angle(V4));
fprintf(tab, "$V_5$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V5),angle(V5));
fprintf(tab, "$V_6$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V6_a),angle(V6_a));
fprintf(tab, "$V_7$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V7),angle(V7));
fprintf(tab, "$V_8$ & $%fe^{%fj}$ \\\\ \\hline \n", abs(V8),angle(V8));
fclose(tab);

%%-----------------> Calculate total solution <-------------------

Gain=abs(V6_a);
Phase=angle(V6_a);

v6_forc=Gain*cos(w*t+Phase);

v6=v6_n+v6_forc;

vs=sin(w*t);

tb=-5e-3:1e-6:0;

V6_b=ones(1,length(tb))*V6_b;
Vs=values(10,3);
Vs=ones(1,length(tb))*Vs;

tot_sol = figure ();
plot (t*1000, v6, "b");  
hold on;
plot (t*1000, vs, "r");
hold on;
plot (tb*1000, V6_b, "b"); 
%hold on;
%plot (0, V6_0, "r"); 
hold on;
plot (tb*1000, Vs, "r");
xlabel ("t [ms]");
ylabel ("Voltage [V]");
legend("v6","vs");
title ("Total Solution");
print (tot_sol, "tot_sol.eps", "-depsc");


%%-----------------> Frequency Response <-------------------

f=logspace(-1,6); %Hz

for i=1:length(f)
	
w=2*pi*f(i);
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

V6(i)=V(5);
V8(i)=V(7);

end

Vs=ones(1,length(f))*Vs;

T6=V6./Vs;
Tc=(V6-V8)./Vs;
Ts=Vs./Vs;
magn6=abs(T6);
magnc=abs(Tc);
magns=abs(Ts);
phase6=angle(T6);
phasec=angle(Tc);
phases=angle(Ts);

freq_db = figure ();
plot (log10(f), 20*log10(magns), "r"); 
hold on;
plot (log10(f), 20*log10(magnc), "y");
hold on;
plot (log10(f), 20*log10(magn6), "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Magnitude v_s(f), v_c(f), v_6(f) [dB]");
legend("v_s(f)","v_c(f)","v_6(f)", 'location', 'southwest');
title ("Frequency Response - Magnitude");
print (freq_db, "freq_db.eps", "-depsc");

phase_ang = figure ();
plot (log10(f), phases*180/pi, "r"); 
hold on;
plot (log10(f), phasec*180/pi, "y");
hold on;
plot (log10(f), phase6*180/pi, "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Phase v_s(f), v_c(f), v_6(f) [degrees]");
legend("v_s(f)","v_c(f)","v_6(f)", 'location', 'southwest');
title ("Frequency Response - Phase");
print (phase_ang, "phase_ang.eps", "-depsc");


%%-----------------> EXPORT TO NGSPICE <-------------------

file=fopen("ngspice_tb0.cir","w");
fprintf(file,".OP\nR1 1 2 %.11fk\nR2 3 2 %.11fk\nR3 2 5 %.11fk\nR4 5 0 %.11fk\nR5 5 6 %.11fk\nR6 0 4 %.11fk\nR7 7 8  %.11fk\nVs 1 0 %.11f\nVo 4 7 0\nC 6 8 %.11fu\nHvd 5 8 Vo %.11fk\nGib 6 3 2 5 %.11fm\n.END\n", values(3,4), values(4,3), values(5,3), values(6,3), values(7,3), values(8,3), values(9,3), values(10,3), values(11,3), values(13,3), values(12,3));
fclose(file);

file=fopen("ngspice_vs0.cir","w");
fprintf(file,".OP\nR1 1 2 %.11fk\nR2 3 2 %.11fk\nR3 2 5 %.11fk\nR4 5 0 %.11fk\nR5 5 6 %.11fk\nR6 0 4 %.11fk\nR7 7 8  %.11fk\nVs 1 0 0\nVo 4 7 0\nVx 6 8 %.11f\nHvd 5 8 Vo %.11fk\nGib 6 3 2 5 %.11fm\n.END\n", values(3,4), values(4,3), values(5,3), values(6,3), values(7,3), values(8,3), values(9,3), Vx, values(13,3), values(12,3));
fclose(file);

file=fopen("ngspice_nat.cir","w");
fprintf(file,".OP\nR1 1 2 %.11fk\nR2 3 2 %.11fk\nR3 2 5 %.11fk\nR4 5 0 %.11fk\nR5 5 6 %.11fk\nR6 0 4 %.11fk\nR7 7 8  %.11fk\nVs 1 0 0\nVo 4 7 0\nC 6 8 %.11fu\nHvd 5 8 Vo %.11fk\nGib 6 3 2 5 %.11fm\n.END\n.ic v(6)=%.11f v(8)=%e\n", values(3,4), values(4,3), values(5,3), values(6,3), values(7,3), values(8,3), values(9,3), values(11,3), values(13,3), values(12,3), V61, V81);
fclose(file);

file=fopen("ngspice_tot.cir","w");
fprintf(file,".OP\nR1 1 2 %.11fk\nR2 3 2 %.11fk\nR3 2 5 %.11fk\nR4 5 0 %.11fk\nR5 5 6 %.11fk\nR6 0 4 %.11fk\nR7 7 8  %.11fk\nVs 1 0 0.0 ac 1.0 sin(0 1 1k)\nVo 4 7 0\nC 6 8 %.11fu\nHvd 5 8 Vo %.11fk\nGib 6 3 2 5 %.11fm\n.ic v(6)=%.11f v(8)=%e\n.END\n", values(3,4), values(4,3), values(5,3), values(6,3), values(7,3), values(8,3), values(9,3), values(11,3), values(13,3), values(12,3), V61, V81);
fclose(file);








