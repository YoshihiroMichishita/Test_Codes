# README

### These are the codes for calculating the nonlinear conductivity with the Green function methods, the reduced density-matrix methods, or the non-Hermitian band-index methods.

　　

　　
　　
---
---


### The main code are made from four object, which are parm_~.o, matrix_op.o, Ham_~.o, and transport.o.

1. parm_~.cpp(hpp) hold the information about the parameter we can input from outside.

1. matrix_op_mypc.cpp(hpp) holds the matrix operation we use in the calculation.

1. Ham_~.cpp(hpp) can set up the Hamiltonian, the velocity operator, and the Green function. In the constructor Ham(parm, k, sw), sw can deside what we use, that is the Green function methods, the reduced density-matrix methods, or the non-Hermitian band-index methods.
When you choose sw=0, you can calculate with the reduced density matix in conventional band-index and not calculate the Green function. When you choose sw=0, you can calculate with the Green function methods. When you choose sw=2, you can calculate with the reduced density matix in non-Hermitian band-index.  

1. transport.cpp(hpp) can calculate the nonlinear conductivity by using Ham~.o.

---


### You can compile the code with makefile.

- CXX: compiler. Usually you can use g++ or icpc.

- LIB: include files. For these codes, you cannot use mkl library.(Because of a bit difference from lapack, and it becomes crutial when you calculate the nonlinear conductivity with band-index.)

- CXX_FLAGS: option for compile. By using O3, the program calculates faster although it takes more time when compiling.

- NAME:the file name for .out file.

- MAIN: the name of main.cpp

- HAM: the name of Ham*.cpp. By changing this, you can change the Hamiltonian you calculate.

- PARM: the name of parm_*.cpp.

- OMP: you choose fopenmp when using g++ for compiler. If you use icpc, you should choose qopenmp. 