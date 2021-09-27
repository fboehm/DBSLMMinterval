# -----------------------------------------------------------------
#   Makefile for dbslmminterval 
# ---------------------------------------------------------------------

# Set the file type
OUTPUTD = bin/dbslmminterval

# Set the path of library
ARMALIB = /net/mulan/home/yasheng/Cpp/arma/lib

# Put C++ complier 
CXX = g++
#CXX = clang++-11

# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): src/main.o src/dtpr.o src/dbslmm.o src/dbslmmfit.o src/calc_asymptotic_variance.o src/subset_to_test_and_training.o
	$(CXX) src/main_dbslmm.o src/dtpr.o src/dbslmm.o src/dbslmmfit.o src/calc_asymptotic_variance.o src/subset_to_test_and_training.o -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)
src/main.o: src/main.cpp
	$(CXX) -c src/main.cpp
src/dbslmm.o: src/dbslmm.cpp include/dbslmm.hpp
	$(CXX) -c src/dbslmm.cpp  $(CXXFLAG)
src/dbslmmfit.o: src/dbslmmfit.cpp include/dbslmmfit.hpp include/calc_asymptotic_variance.h include/subset_to_test_and_training.h
	$(CXX) -c src/dbslmmfit.cpp $(CXXFLAG)
src/dtpr.o: src/dtpr.cpp include/dtpr.hpp 
	$(CXX) -c src/dtpr.cpp $(CXXFLAG)

src/calc_asymptotic_variance.o: src/calc_asymptotic_variance.cpp include/calc_asymptotic_variance.h 
	$(CXX) -c src/calc_asymptotic_variance.cpp $(CXXFLAG)
src/subset_to_test_and_training.o: src/subset_to_test_and_training.cpp include/subset_to_test_and_training.h 
	$(CXX) -c src/subset_to_test_and_training.cpp $(CXXFLAG)



clean:
	rm -f src/*.o  bin/dbslmm
