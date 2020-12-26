SRC=$(PWD)/src
OBJ=$(PWD)/obj

CXX=icpc
#INC=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
LIB=-llapack -lblas
#LIB=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#CXX_FLAGS=-O3 -I$(INC) -L$(LIB)
CXX_FLAGS=-O3 $(LIB)

NAME=TMD_BI
MAIN=main_BI
HAM=Ham_TMD
MATRIX=matrix_op_mypc

### make all
all: $(NAME)

$(NAME): $(OBJ)/$(MAIN).o $(OBJ)/C3v_HSL.o $(OBJ)/transport.o $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/parm_TMD.o
	$(CXX) -qopenmp $(OBJ)/$(MAIN).o $(OBJ)/C3v_HSL.o $(OBJ)/transport.o $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/parm_TMD.o $(CXX_FLAGS) -o $(NAME) 

$(OBJ)/$(MAIN).o: $(SRC)/$(MAIN).cpp
	$(CXX) -qopenmp -c $(SRC)/$(MAIN).cpp -o $(OBJ)/$(MAIN).o $(CXX_FLAGS)

$(OBJ)/C3v_HSL.o: $(SRC)/C3v_HSL.cpp $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/parm_TMD.o
	$(CXX) -c $(SRC)/C3v_HSL.cpp -o $(OBJ)/C3v_HSL.o  $(CXX_FLAGS)

$(OBJ)/transport.o: $(SRC)/transport.cpp $(OBJ)/$(HAM).o $(OBJ)/$(MATRIX).o $(OBJ)/parm_TMD.o
	$(CXX) -c $(SRC)/transport.cpp -o $(OBJ)/transport.o  $(CXX_FLAGS)

$(OBJ)/$(HAM).o: $(SRC)/$(HAM).cpp $(OBJ)/$(MATRIX).o $(OBJ)/parm_TMD.o
	$(CXX) -c $(SRC)/$(HAM).cpp -o $(OBJ)/$(HAM).o $(CXX_FLAGS)

$(OBJ)/$(MATRIX).o: $(SRC)/$(MATRIX).cpp
	$(CXX) -c $(SRC)/$(MATRIX).cpp -o $(OBJ)/$(MATRIX).o $(CXX_FLAGS)

$(OBJ)/parm_TMD.o: $(SRC)/parm_TMD.cpp
	$(CXX) -c $(SRC)/parm_TMD.cpp -o $(OBJ)/parm_TMD.o $(CXX_FLAGS)

clean:
	rm $(OBJ)/* $(PWD)/a.out

