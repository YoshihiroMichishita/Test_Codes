SRC=$(PWD)/src
OBJ=$(PWD)/obj

CXX=g++
#INC=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
LIB=-llapack -lblas -lgomp -lpthread -lm -ldl
#LIB=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#CXX_FLAGS=-O3 -I$(INC) -L$(LIB)
CXX_FLAGS=-O3 $(LIB)

NAME=Dirac_BI_Green
MAIN=main_DiracG2_mypc
HAM=Ham_Dirac
MATRIX=matrix_op_mypc
PARM=parm_Dirac
OMP=fopenmp

### make all
all: $(NAME)

$(NAME): $(OBJ)/$(MAIN).o $(OBJ)/transport.o $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/$(PARM).o
	$(CXX) -$(OMP) $(OBJ)/$(MAIN).o $(OBJ)/transport.o $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/$(PARM).o $(CXX_FLAGS) -o $(NAME) 

$(OBJ)/$(MAIN).o: $(SRC)/$(MAIN).cpp
	$(CXX) -$(OMP) -c $(SRC)/$(MAIN).cpp -o $(OBJ)/$(MAIN).o $(CXX_FLAGS)

#$(OBJ)/C3v_HSL.o: $(SRC)/C3v_HSL.cpp $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/$(PARM).o
#	$(CXX) -c $(SRC)/C3v_HSL.cpp -o $(OBJ)/C3v_HSL.o  $(CXX_FLAGS)

$(OBJ)/transport.o: $(SRC)/transport.cpp $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/$(PARM).o
	$(CXX) -c $(SRC)/transport.cpp -o $(OBJ)/transport.o  $(CXX_FLAGS)

$(OBJ)/$(HAM).o: $(SRC)/$(HAM).cpp $(OBJ)/$(MATRIX).o $(OBJ)/$(PARM).o
	$(CXX) -c $(SRC)/$(HAM).cpp -o $(OBJ)/$(HAM).o $(CXX_FLAGS)

$(OBJ)/$(MATRIX).o: $(SRC)/$(MATRIX).cpp
	$(CXX) -c $(SRC)/$(MATRIX).cpp -o $(OBJ)/$(MATRIX).o $(CXX_FLAGS)

$(OBJ)/$(PARM).o: $(SRC)/$(PARM).cpp
	$(CXX) -c $(SRC)/$(PARM).cpp -o $(OBJ)/$(PARM).o $(CXX_FLAGS)

clean:
	rm $(OBJ)/* $(PWD)/a.out

