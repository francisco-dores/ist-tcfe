# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

number=95790

all:	python
	make -C mat
	make -C sim
	make -C doc
	cp doc/report.pdf .

python:
	echo "$(number)" | python2.7 t2_datagen.py > data.txt

clean:
	make -C mat clean
	make -C sim clean
	make -C doc clean
	@rm -rf *.txt

.PHONY: all clean python

