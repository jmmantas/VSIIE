LISROOT     = usr/local
SRC         = src
INC		 = include
INC_Solvers         = ${INC}/Solvers
INC_IVP     = ${INC}/IVPs
DIR_VSIIE    = ${SRC}/VSIIE
DIR_VSBDF    = ${SRC}/VSBDF
OBJ_SOLVER   = IVP_Solver
CHECK_ORDER  = check_order
VSIIE_Solver = VSIIE_Solver
SBDF_Solver = SBDF_Solver
UTILS=utils

IVP=IVP_ODE
IVP1=IVP_ODE_advdiff1d
IVP2=IVP_ODE_stiff_brusselator
IVP3=IVP_ODE_combustion
IVP4=IVP_ODE_HIRES

ALL_IVPs=${INC_IVP}/$(IVP1).h ${INC_IVP}/$(IVP2).h ${INC_IVP}/$(IVP3).h ${INC_IVP}/$(IVP4).h   ${INC_IVP}/$(IVP).h

CCLINKLIBS = -L$(LISROOT)/lib  -llis -lopenblas

INCLUDE= -I$(LISROOT)/include/ -I${INC} -I${INC_IVP} -I${INC_Solvers}

# C++ Compiler:
CC=g++

# General compiler flags:
C_FLAGS= -I. ${INCLUDE} -std=c++11  -Wall  -O3 -m64 -ffast-math -fomit-frame-pointer -fpermissive -fopenmp -g

all: ${OBJ_SOLVER}  ${CHECK_ORDER}

${OBJ_SOLVER}: ${SRC}/$(OBJ_SOLVER).cpp ${INC}/${UTILS}.h ${UTILS}.o  ${ALL_IVPs}   ${VSIIE_Solver}.o ${SBDF_Solver}.o  
	$(CC)   $(C_FLAGS) ${SRC}/$(OBJ_SOLVER).cpp ${VSIIE_Solver}.o ${SBDF_Solver}.o ${UTILS}.o $(CCLINKLIBS) -o $(OBJ_SOLVER) 

${CHECK_ORDER}: ${SRC}/$(CHECK_ORDER).cpp ${INC}/${UTILS}.h  ${UTILS}.o ${ALL_IVPs}   ${VSIIE_Solver}.o   
	$(CC)  $(C_FLAGS)   ${SRC}/$(CHECK_ORDER).cpp ${VSIIE_Solver}.o ${UTILS}.o $(CCLINKLIBS) -o $(CHECK_ORDER) 

${VSIIE_Solver}.o:  ${DIR_VSIIE}/${VSIIE_Solver}.cpp ${INC}/${UTILS}.h  ${UTILS}.o 
	$(CC) -c  $(C_FLAGS)  ${DIR_VSIIE}/${VSIIE_Solver}.cpp  -o  ${VSIIE_Solver}.o 

${SBDF_Solver}.o:  ${DIR_VSBDF}/${SBDF_Solver}.cpp  ${INC}/${UTILS}.h ${UTILS}.o 
	$(CC) -c  $(C_FLAGS) ${DIR_VSBDF}/${SBDF_Solver}.cpp  -o  ${SBDF_Solver}.o 
	
${UTILS}.o:  ${INC}/${UTILS}.h   
	$(CC) -c  $(C_FLAGS) ${SRC}/$(UTILS).cpp  -o  ${UTILS}.o	




.PHONY: clean

clean:
	rm  $(OBJ_SOLVER)  *.o 
cleanout:
	rm  *.txt 
