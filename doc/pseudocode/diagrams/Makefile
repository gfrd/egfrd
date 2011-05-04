
# sudo yum install graphviz

all: egfrd fireSingle firePair fireMulti createMulti

egfrd:
	dot egfrd.dot -Tpng -o egfrd.png

fireSingle:
	dot fireSingle.dot -Tpng -o fireSingle.png

firePair:
	dot firePair.dot -Tpng -o firePair.png

fireMulti:
	dot fireMulti.dot -Tpng -o fireMulti.png

createMulti:
	dot createMulti.dot -Tpng -o createMulti.png

clean:
	rm -rf *.png
