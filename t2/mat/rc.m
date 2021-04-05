close all
clear all


vs(t) = 0 %no excitation

wn=-1/R*C


hf = figure ();
plot (t*1000, vs, "g");
hold on;
plot (t*1000, vcn, "b");

xlabel ("t[ms]");
ylabel ("vs(t), vcn(t) [V]");
print (hf, "forced.eps", "-depsc");
