%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
RE1=150
RC1=2400
RB1=79000
RB2=9200
VBEON=0.7
VCC=12
RS=100

RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1

VCE>VBEON

gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

%AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
%AVI_DB = 20*log10(abs(AV1))
%AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
%AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
%AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
%AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=150
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
%ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

%RE1=0
%ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
%ZO1 = 1/(1/ro1+1/RC1)


tab=fopen("op1.tex", "w");
fprintf(tab, "$v_{CE}$ & $%f$ \\\\ \\hline \n", 20*log10(abs(VCE)));
fprintf(tab, "$v_{O1}$ & $%f$ \\\\ \\hline \n", 20*log10(abs(VO1)));
fprintf(tab, "$I_{C1}$ & $%f$ \\\\ \\hline \n", IC1);
fprintf(tab, "$I_{E1}$ & $%f$ \\\\ \\hline \n", IE1);
fprintf(tab, "$I_{B1}$ & $%f$ \\\\ \\hline \n", IB1);
fclose(tab);

tab=fopen("imp_gain1.tex", "w");
fprintf(tab, "$Z_{I1}$ & $%f$ \\\\ \\hline \n", ZI1);
fprintf(tab, "$Z_{O1}$ & $%f$ \\\\ \\hline \n", ZO1);
fprintf(tab, "$A_{V1}$ & $%f$ \\\\ \\hline \n", 20*log10(abs(AV1)));
fclose(tab);



%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 320
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)


tab=fopen("op2.tex", "w");
fprintf(tab, "$v_{O2}$ & $%f$ \\\\ \\hline \n", 20*log10(abs(VO2)));
fprintf(tab, "$I_{C2}$ & $%f$ \\\\ \\hline \n", IC2);
fprintf(tab, "$I_{E2}$ & $%f$ \\\\ \\hline \n", IE2);
fprintf(tab, "$I_{B2}$ & $%f$ \\\\ \\hline \n", IB1);
fclose(tab);

tab=fopen("imp_gain2.tex", "w");
fprintf(tab, "$Z_{I2}$ & $%f$ \\\\ \\hline \n", ZI2);
fprintf(tab, "$Z_{O2}$ & $%f$ \\\\ \\hline \n", ZO2);
fprintf(tab, "$A_{V2}$ & $%f$ \\\\ \\hline \n", 20*log10(abs(AV2)));
fclose(tab);

tab=fopen("imp_gaintot.tex", "w");
fprintf(tab, "$Z_I$ & $%f$ \\\\ \\hline \n", ZI);
fprintf(tab, "$Z_O$ & $%f$ \\\\ \\hline \n", ZO);
fprintf(tab, "$A_V$ & $%f$ \\\\ \\hline \n", AV_DB);
fclose(tab);
