SHELL:=/bin/bash
# Makefile

MERGER=file-merger

DECODER=thrsana
DEPS=utils.h progressbar.h

#targets
all: $(MERGER) $(DECODER)

$(MERGER): $(MERGER).cpp
	g++ --std=c++17 -march=native -gdwarf-4 -pedantic -Wall -O3 -o $(MERGER) $(MERGER).cpp -lm -g

$(DECODER): $(DECODER).cpp $(DEPS)
	g++ --std=c++17 -march=native -mavx2 -gdwarf-4 -pedantic -Wall -O3 -o $(DECODER) $(DECODER).cpp -lm -g


clean:
	rm -f $(MERGER) $(DECODER)

help:
	(cat Makefile | grep '##')

%:
	@:
