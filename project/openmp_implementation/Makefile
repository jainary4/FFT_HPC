# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++14 -O2 -Wall 

# OpenMP flag
OPENMP_FLAG = -fopenmp

# Profiling flags
PROFILING_FLAGS = -pg

# Target executable
TARGET = fft

# Source files
SRC = fft.cpp \
      complex.h \
      dit_fft.h \
      inverse_dit.h \
      denoise.h \

# Object files (only fft.cpp needs compilation)
OBJ = fft.o

# Default target (no OpenMP)
all: $(TARGET)

# Build the target without OpenMP
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Compile fft.cpp into an object file
fft.o: fft.cpp complex.h dit_fft.h inverse_dit.h validate_n_evaluate.h denoise.h derivative.h
	$(CXX) $(CXXFLAGS) -c fft.cpp -o fft.o

# Parallel target (with OpenMP)
parallel: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) $(OPENMP_FLAG)" $(TARGET)

# Profiling target for building with profiling flags and OpenMP
profile: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) $(OPENMP_FLAG) $(PROFILING_FLAGS)" $(TARGET)

# Clean up build files
clean:
	rm -f $(OBJ) $(TARGET) gmon.out
	rm -f fft_output.txt
	rm -f ifft_output.txt

# Run the program
run: $(TARGET)
	 ./$(TARGET) 25  3

