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
fprintf(tab, "$f_L$ [Hz] & $%f$ \\\\ \\hline \n", f_low);
fprintf(tab, "$f_H$ [Hz] & $%f$ \\\\ \\hline \n", f_high);
fprintf(tab, "$f_0$ [Hz] & $%f$ \\\\ \\hline \n", f_c);
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
phase_Av_deg = angle(T)*180/pi;
%phase_Av_rad = angle(T);

end


freq_db = figure ();
plot (log10(f), 20*log10(Av), "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Magnitude [dB]");
title ("Gain - Magnitude");
print (freq_db, "freq_db.eps", "-depsc");

freq_p_deg = figure ();
plot (log10(f), phase_Av_deg, "r");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Phase [deg]");
title ("Gain - Phase");
print (freq_p_deg, "freq_p_deg.eps", "-depsc");

tab=fopen("imped.tex", "w");
fprintf(tab, "$Z_I$ [S] & $%f$ \\\\ \\hline \n", ZI);
fprintf(tab, "$Z_O$ [S] & $%f$ \\\\ \\hline \n", ZO);
fclose(tab);

tab=fopen("freq_dev.tex", "w");
fprintf(tab, "$f_0 deviation$ [Hz] & $%f$ \\\\ \\hline \n", freq_dev);
fclose(tab);

tab=fopen("gain.tex", "w");
fprintf(tab, "$A_v$ [dB] & $%f$ \\\\ \\hline \n", Av_1_db);
fprintf(tab, "$A_v deviation$ [dB] & $%f$ \\\\ \\hline \n", gain_dev_db);
fclose(tab);

