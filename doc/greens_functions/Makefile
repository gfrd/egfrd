
.PHONY: all
all:
	cd tex && $(MAKE) $(MKFLAGS) && cd ../ && mv -f build/*.pdf .

clean:
	rm -rf build/*.aux build/*.log build/include/*.aux build/include/*.log *.pdf
