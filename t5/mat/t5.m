R1 = 1000;
R2 = 1e3;
R3 = 300e3;
R4 = 5e3;
C1 = 220e-9;
C2 = 110e-9;


f_low = 1/(2*pi*C1*R1);
f_high = 1/(2*pi*C2*R2);

%CENTRAL FREQUENCY

f_c=sqrt(f_low*f_high)

tab=fopen("freq.tex", "w");
fprintf(tab, "$f_L$ & $%f$ \\\\ \\hline \n", f_low);
fprintf(tab, "$f_H$ & $%f$ \\\\ \\hline \n", f_high);
fprintf(tab, "$f_0$ & $%f$ \\\\ \\hline \n", f_c);
fclose(tab);


%CENTRAL FREQUENCY DEVIATION
cfreq_dev = abs(f_c-1000)


%TRANSFER FUNCTION

%T = ((R1*C1*s)/(1+R1*C1*s))*(1+R3/R4)*(1/(1+R2*C2*s));


%INPUT IMPEDANCE (f=f_c)

ZI = R1 + 1/(2*pi*f_c*C1)


%OUTPUT IMPEDANCE (f=f_c)

ZO = 1/(2*pi*f_c*C2+1/R2)


%GAIN (f=f_c)

w=2*pi*f_c;
s=j*w;
T = ((R1*C1*s)/(1+R1*C1*s))*(1+R3/R4)*(1/(1+R2*C2*s));

Av_c = abs(T)


%GAIN (f=1kHz)

w=2*pi*1000;
s=j*w;
T = ((R1*C1*s)/(1+R1*C1*s))*(1+R3/R4)*(1/(1+R2*C2*s));

Av_1 = abs(T)


%GAIN DEVIATION
gain_dev = abs(Av_1-100)
 
 
%GAIN DEVIATION dB
Av_1_db = 20*log10(Av_1)
gain_dev_db = abs(Av_1_db-40)

%FREQUENCY DEVIATION Hz

freq_dev=1000-f_c;


%FREQUENCY RESPONSE

f=logspace(1,8); %Hz
for i=1:length(f)
w=2*pi*f(i);
s = j*w;

T(i) = ((R1*C1*s)/(1+R1*C1*s))*(1+R3/R4)*(1/(1+R2*C2*s));
Av = abs(T);
phase_Av = angle(T)*180/pi;

end


freq_db = figure ();
plot (log10(f), 20*log10(Av), "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Gain - Magnitude [dB]");
title ("Gain - Magnitude");
print (freq_db, "freq_db.eps", "-depsc");

freq_p = figure ();
plot (log10(f), phase_Av, "r");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Gain - Phase [deg]");
title ("Gain - Phase");
print (freq_p, "freq_p.eps", "-depsc");


tab=fopen("imped.tex", "w");
fprintf(tab, "$Z_I$ & $%f$ \\\\ \\hline \n", ZI);
fprintf(tab, "$Z_O$ & $%f$ \\\\ \\hline \n", ZO);
fclose(tab);

tab=fopen("dev.tex", "w");
fprintf(tab, "$A_v deviation (Db)$ & $%f$ \\\\ \\hline \n", gain_dev_db);
fprintf(tab, "$f_0 deviation (Hz)$ & $%f$ \\\\ \\hline \n", freq_dev);
fclose(tab);



