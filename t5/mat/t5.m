R1 = 500;
R2 = 1e3;
R3 = 300e3;
R4 = 5e3;
C1 = 220e-9;
C2 = 110e-9;


f_low = 1/(2*pi*C1*R1);
f_high = 1/(2*pi*C2*R2);

%CENTRAL FREQUENCY

f_c=sqrt(f_low*f_high)


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


%FREQUENCY RESPONSE

f=logspace(1,8); %Hz
for i=1:length(f)
w=2*pi*f(i);
s = j*w;

T(i) = ((R1*C1*s)/(1+R1*C1*s))*(1+R3/R4)*(1/(1+R2*C2*s));
Av = abs(T);


end


freq_db = figure ();
plot (log10(f), 20*log10(Av), "b");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Gain [dB]");
title ("Frequency Response - Gain");
print (freq_db, "freq_db.eps", "-depsc");



