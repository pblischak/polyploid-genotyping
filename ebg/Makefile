# Makefile for compiling
EXE = ebg
OBJ = MbRandom.o ModelGeneric.o \
			ModelHWE.o ModelDiseq.o \
			ModelAlloSNP.o ModelAlloSNP2.o \
			ModelGATK.o ModelGL.o main.o
CXX = g++
CXXFLAGS = -Wall -O3 -g -std=c++11 -Wextra -Wpedantic
OPENMP=no
ifeq ($(strip $(OPENMP)), yes)
        CXXFLAGS += -fopenmp
endif
.PHONY : clean install uninstall

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
	@printf "\nTo install ebg system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files and executable...\n\n"
	@rm $(OBJ) $(EXE)

install :
	@printf "\nCopying executable to /usr/local/bin...\n\n"
	@cp $(EXE) /usr/local/bin
	@printf "\nTo uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\nRemoving executable from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE)

test :
	./ebg hwe -n 100 -l 10 -p 4 -t tot.txt -a alt.txt -e err.txt
