CXXFLAGS := -O2 -g -Wall -std=c++11
BINARY := branchsim

all: $(BINARY)

branchsim: $(BINARY).o branchsim_driver.o
	$(CXX) -o $(BINARY) -lm branchsim.o branchsim_driver.o


clean:
	rm -rf $(BINARY) *.o myoutput
