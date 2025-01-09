# Compiler and flags
CXX = mpicxx
CXXFLAGS = -std=c++14 -O2 -Wall

# OpenMP flag
OPENMP_FLAG = -fopenmp

# Profiling flags
PROFILING_FLAGS = -pg

# Target executable
TARGET = fft_hybrid

# Source files
SRC = fft.cpp \
      complex.h \
      dit_fft.h \
      inverse_dit.h \
      denoise.h

# Object files (only fft.cpp needs compilation)
OBJ = fft.o

# Default target (MPI + OpenMP)
all: $(TARGET)

# Build the hybrid target (MPI + OpenMP)
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) -o $(TARGET) $(OBJ)

# Compile fft.cpp into an object file
fft.o: fft.cpp complex.h dit_fft.h inverse_dit.h validate_n_evaluate.h denoise.h derivative.h
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) -c fft.cpp -o fft.o

# Profiling target for building with profiling flags and hybrid MPI + OpenMP
profile: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) $(OPENMP_FLAG) $(PROFILING_FLAGS)" $(TARGET)

# Clean up build files
clean:
	rm -f $(OBJ) $(TARGET) gmon.out
	rm -f fft_output.txt
	rm -f ifft_output.txt

# Run the program with MPI + OpenMP
run: $(TARGET)
	mpirun -np 16 ./$(TARGET) 25 3
	

