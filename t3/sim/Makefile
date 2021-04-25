all: t2_tb0 t2_vs0 t2_natsolution t2_totalsolution

t2_tb0: t2_tb0.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $ttb0_tab.tex;)

t2_tb0.log: ../sim/t2_tb0.net
	ngspice -b $< -o $@

t2_vs0: t2_vs0.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $vvs0_tab.tex;)

t2_vs0.log: ../sim/t2_vs0.net
	ngspice -b $< -o $@

t2_natsolution: t2_natsolution.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)

t2_natsolution.log: ../sim/t2_natsolution.net
	ngspice -b $< -o $@

t2_totalsolution: t2_totalsolution.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)

t2_totalsolution.log: ../sim/t2_totalsolution.net
	ngspice -b $< -o $@

clean:
	@rm -f t2_tb0.log t2_vs0.log t2_natsolution.log t2_totalsolution.log *.tex *.pdf *.ps *~

.PHONY: all clean
