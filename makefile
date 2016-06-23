VRNA := /home/kale/research/software/third_party/ViennaRNA-2.2.5

CC := g++
SRCDIR := src
BUILDDIR := build
TARGET := bin/mh
 
SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -Werror -std=c++11 -fopenmp
LIB := -L/home/kale/.local/lib -lRNA -fopenmp
INC := -isystem /home/kale/.local/include -isystem /usr/include/eigen3

run: $(TARGET)
	$(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $^ -o $(TARGET) $(LIB) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	$(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean run
