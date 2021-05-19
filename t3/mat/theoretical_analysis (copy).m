close all 
clear all


var = true;
num_iteracoes=0; % algorithm iterations
counter=0;
while(var == true)
    f = [];
    alfa = []; % intensidade do passo
    
    
    f(1) = calcula_merito(R1,C,R2,AA);
    
    grad_f =  calcula_grad(R1,C,R2,AA);
    s_ini = grad_f;
    
    s = s_ini/(sqrt(s_ini(1)^2+s_ini(2)^2+s_ini(3)^2+s_ini(4)^2)); % s normalizado
    
​
    var2 = true;
    alfa(1) = 0.05;
    j = 1;
    f_comp = f(1); % f para comparar
    passo_bom = 0;
    while (var2 == true)    
             
        f(j+1) = calcula_merito(R1,C,R2,AA);
        
        if (f(j+1) < f_comp)
            alfa(j+1) = alfa(j)*2;
            f_comp = f(j+1);
            passo_bom = 1;
        elseif (f(j+1) > f_comp && j > 3 && passo_bom == 1)
           
            [~,ind] = min(f);
​
            
            alpha_opt = alfa(ind-1);
            x = x + alpha_opt*s;   
            var2 = false;
        elseif 1e-6 > alfa(j)   
            var2 = false;
            
        else
            alfa(j+1) = alfa(j)*0.8;
        end
        j = j+1;
        x_matrix(counter,:)=[x_new,f_comp];
    end
    
   
    f_new = calcula_merito(R1,C,R2,AA);
    if(abs(f_new - f(1)) < 1e-4)
        var = false;
    end 
    num_iteracoes=num_iteracoes+1;
   
    x_matrix(counter,:)=[x,f_new];
    
end
​
function [grad_f] = calcula_grad(R1,C,R2,AA)

    % primeiro ponto
    [f_mais]= calcula_merito(R1+delta(1),C,R2,AA);
    [f_menos]=calcula_merito(R1-delta(1),C,R2,AA);
    grad_f(1) = (f_mais-f_menos)/(2*delta(1)); 
    
    % segundo ponto
    
     [f_mais]= calcula_merito(R1,C+delta(2),R2,AA);
    [f_menos]=calcula_merito(R1,C-delta(2),C,R2,AA);
    grad_f(2) = (f_mais-f_menos)/(2*delta(2)); 
    
    % terceiro ponto
     [f_mais]= calcula_merito(R1,C,R2+delta(3),AA);
    [f_menos]=calcula_merito(R1,C,R2-delta(3),AA);
    grad_f(3) = (f_mais-f_menos)/(2*delta(3)); 
    
     % quarto ponto
     [f_mais]= calcula_merito(R1,C,R2,AA+delta(4));
    [f_menos]=calcula_merito(R1,C,R2,AA-delta(4));
    grad_f(4) = (f_mais-f_menos)/(2*delta(4)); 
    
end










function [M]=calcula_merito(R1,C,R2,AA)
f=50;
w=2*pi*f;
T=1/f;

t=linspace(0, 10*T/2, 4000);

AS=230;
vS=AS*cos(w*t);

%------------------Transformer-------------------------
%AA=50;
vA=AA*cos(w*t);
n_transformer=AS/AA;

%------------------Full wave rectifier-----------------
vB = abs(vA);

%------------------Envelope detector------------------- 
%R1=1000;
%C=0.1e-3;

TB=T/2;
wB=w*2;

%Diode OFF
tOFF=1/wB*atan(1/(wB*R1*C));


%v_ripple=AA*(1-exp(-TB/(R1*C)));

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

ripple_env=max(vC)-min(vC);

plot(t*1000, vB)
hold on
plot(t*1000,vC)
title("Envelpe detector")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc")
hold off


%--------------Voltage regulator----------------------

%Incremental analysis
%vC = VC + vc
%VON=0.7;
%R2=10000;
%num_diodes=round(12/VON);
num_diodes=17;
IS=1e-14;
VT=25e-3;
eta=1;
VC= sum(vC)/length(vC);

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
%disp(VO(end));
VON=VO(end)/num_diodes

rd=eta*VT/(IS*exp(VON/(eta*VT))); %ver bem isto


for i=1:length(vC)
vo(i)=(num_diodes*rd/(R2+num_diodes*rd))*(vC(i)-VC);
end


VO=num_diodes*VON; %valor do método Newton-Rhapson
vO=VO+vo;

plot(t*1000, vO)
hold on
plot(t*1000,vC)
title("")
xlabel ("t[ms]")
legend("output voltage","envelope")
print ("vregulator.eps", "-depsc")
hold off

plot(t*1000,vO-12)
title("Deviation")
xlabel ("t[ms]")
legend("vO-12")
print ("vdeviation.eps", "-depsc")


ripple_out= max(vO)-min(vO)
DC_out=sum(vO)/length(vO)

M=1/(((R1+R2)/1000+C/1e6)*(ripple_out+abs(DC_out-12)+10e-6));
end





