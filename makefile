CXX := g++
SRCEXT := cc
CFLAGS := -g -Wall -Werror -std=c++11 -fopenmp
LIB := -L/home/kale/.local/lib \
       -L/usr/local/lib64 \
       -lRNA \
       -l:libdocopt.a
INC := -isystem /home/kale/.local/include \
       -isystem /usr/local/include \
       -I include

SRC_OBJS := $(patsubst %.cc,build/%.o,$(wildcard src/*.cc))
TEST_OBJS := $(patsubst %.cc,build/%.o,$(wildcard tests/*.cc))

# Rules for building apps.

mh: bin/mh
	$< -n100 2> bin/mh.log

bin/%: build/apps/%.o $(SRC_OBJS) $(EXTERNAL_OBJS)
	@mkdir -p bin
	$(CXX) $(CFLAGS) $^ -o $@ $(LIB)

build/apps/%.o: apps/%.cc
	@mkdir -p build/apps
	$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

build/src/%.o: src/%.cc include/sgrna_design/%.hh
	@mkdir -p build/src
	$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

# Rules for building tests.

test: tests/run_tests
	$<

tests/run_tests: $(TEST_OBJS) $(SRC_OBJS)
	$(CXX) $(CFLAGS) $^ -o $@ $(LIB)

build/tests/%.o: tests/%.cc
	@mkdir -p build/tests
	$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

# Clean up.

clean:
	$(RM) -r build bin tests/run_tests

.PHONY: mh test clean
