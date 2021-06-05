RI=1000;
RF=1000;
CI=220e-9;
CF=220e-9;


f_low = 1/(2*pi*CF*RF);
f_high = 1/(2*pi*CI*RI);

f_c=sqrt(f_low*f_high);


%Input Impedance (f=f_c)

ZI = RI + 1/(2*pi*f_c*CI);


%Output Impedance (f=f_c)

ZF = (RF*1/(2*pi*CF))/(RF+1/(2*pi*CF));


%Gain (f=f_c)

Av = -ZF/ZI;



%Frequency Response

f=logspace(1,8); %Hz
for i=1:length(f)
w=2*pi*f(i);

ZI = RI + 1/(j*w*CI);
ZF = (RF*1/(j*w*CF))/(RF+1/(j*w*CF)); %QUAL DAS DUAS FÓRMULAS??? ESTA OU...
Av (i) = -ZF/ZI; 

end


freq_db = figure ();
plot (log10(f), 20*log10(Av), "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Gain [dB]");
title ("Frequency Response - Gain");
print (freq_db, "freq_db.eps", "-depsc");



