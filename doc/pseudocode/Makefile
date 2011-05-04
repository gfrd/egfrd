
# sudo yum install python-pygments

all: pseudocode diagrams

pseudocode:
	pygmentize -l python -f html -O full -o pseudocode.html pseudocode.py

.PHONY: diagrams
diagrams:
	cd diagrams && make

clean:
	rm -rf pseudocode.html
	cd diagrams && make clean
