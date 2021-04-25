close all 
clear all

f=50;
w=2*pi*f;
T=1/f;

t=linspace(0, 10*T/2, 4000);

AS=230;
vS=AS*cos(w*t);

%------------------Transformer-------------------------
AA=12;
vA=AA*cos(w*t);
n_transformer=AS/AA;

%------------------Full wave rectifier-----------------
vB = abs(vA);

%------------------Envelope detector------------------- 
R1=1000;
C=1e-3;

TB=T/2;
wB=w*2;

%Diode OFF
tOFF=1/wB*atan(1/(wB*R1*C));


v_ripple=AA*(1-exp(-TB/(R1*C)));
j=1;
for ciclo=1:10
	for i=1:length(t)
		%Diode ON
		if t(i)>=(ciclo-1)*TB && t(i)<=(ciclo-1)*TB+tOFF
			vC(i)=vB(i);
		end
		
		%Diode OFF
		v_exp=AA*cos(wB*((ciclo-1)*TB+tOFF))*exp(-(t(i)-((ciclo-1)*TB+tOFF))/(R1*C));
		if t(i)>((ciclo-1)*TB+tOFF) && v_exp>=vB(i) && t(i)<=(ciclo)*TB
			vC(i)=v_exp;
		end
		
		%Diode ON
		if vB(i)>v_exp && t(i)<=(ciclo)*TB
			vC(i)=vB(i);
		end	
	end
end

plot(t*1000, vB)
hold on
plot(t*1000,vC)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc")


%--------------Voltage regulator----------------------

%Incremental analysis
%vC = VC + vc
VON=0.7;
R2=3000;

IS=1e-14;
VT=25e-3;
eta=1;

VC= sum(vC)/length(vC);
rd=eta*VT/IS/exp(VC/eta/VT);

num_diodes=round(12/VON);

VO=num_diodes*VON;

vo=(num_diodes*rd/(R2+num_diodes*rd))*v_ripple;

f=@(VO) -VC+IS*exp(VO/num_diodes/eta/VT)*R2+VO;
h=0.001;
i=1;
VO(i)=12;
err=1;
while err>0.0001

	dif_finita=(f(VO(i)+h)-f(VO(i)-h))/2/h;
	VO(i+1)=VO(i)-f(VO(i))/dif_finita;
	err=VO(i+1)-VO(i);
	i=i+1;
end

disp(VO(end));










