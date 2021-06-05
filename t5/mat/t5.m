RI=;
RF=;
CI=;
CF=;


f_low = 1/(2*pi*CF*RF);
f_high = 1/(2*pi*CI*RI);

f_c=sqrt(f1*f2);


%Input Impedance (f=f_c)

ZI = RI + 1/(2*pi*f_c*CI);


%Output Impedance (f=f_c)

ZF = 1/(1/RF+1/1/(2*pi*f_c*CF));


%Gain (f=f_c)

A_v = -ZF/ZI



%Frequency Response


