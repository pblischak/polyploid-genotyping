# Makefile for compiling

EXE = ebg
OBJ = MbRandom.o ModelGeneric.o ModelHWE.o ModelDiseq.o ModelAlloSNP.o main.o

CXX = g++
CXXFLAGS = -Wall -O3 -g -std=c++11 -march=native

.PHONY : clean install uninstall

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
	@printf "\nTo install ebg system-wide run: sudo make install\n\n"

#linux : $(OBJ)
#	$(CXX) $(CXXFLAGS) -L${BOOST_LIB} -Wl,-rpath=${BOOST_LIB} -o $(EXE) $^ $(LDFLAGS)
#	@printf "\nTo install ppgtk system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files and executable...\n\n"
	@rm -i $(OBJ) $(TXT_FILES) $(EXE)

install :
	@printf "\nCopying executable to /usr/local/bin...\n\n"
	@cp $(EXE) /usr/local/bin
	@printf "\nTo uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\nRemoving executable from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE)
