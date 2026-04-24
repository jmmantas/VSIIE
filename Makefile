LISROOT     = usr/local
SRC         = src
INC         = include
DIR_IVP     = ${INC}/IVPs
DIR_VSIIE    = ${SRC}/VSIIE
OBJ_SOLVER   = IVP_Solver
CHECK_ORDER  = check_order
VSIIE_Solver = VSIIE_Solver


IVP=IVP_ODE
IVP1=IVP_ODE_advdiff1d
IVP2=IVP_ODE_stiff_brusselator
IVP3=IVP_ODE_combustion
IVP4=IVP_ODE_HIRES
ALL_IVPs=${DIR_IVP}/$(IVP1).h ${DIR_IVP}/$(IVP2).h ${DIR_IVP}/$(IVP3).h ${DIR_IVP}/$(IVP4).h   ${DIR_IVP}/$(IVP).h

CCLINKLIBS = -L$(LISROOT)/lib  -llis -lopenblas

# C++ Compiler:
CC=g++

# General compiler flags:
C_FLAGS= -I. -I$(LISROOT)/include/ -std=c++11  -Wall  -O3 -m64 -ffast-math -fomit-frame-pointer -fpermissive -fopenmp -g

all: ${OBJ_SOLVER}  ${CHECK_ORDER}

${OBJ_SOLVER}: ${SRC}/$(OBJ_SOLVER).cpp ${INC}/utils.h  ${ALL_IVPs}   ${VSIIE_Solver}.o   
	$(CC)   $(C_FLAGS) -I${DIR_IVP} -I${INC} ${SRC}/$(OBJ_SOLVER).cpp ${VSIIE_Solver}.o  $(CCLINKLIBS) -o $(OBJ_SOLVER) 

${CHECK_ORDER}: ${SRC}/$(CHECK_ORDER).cpp ${INC}/utils.h  ${ALL_IVPs}   ${VSIIE_Solver}.o   
	$(CC)  $(C_FLAGS)  -I${DIR_IVP} -I${INC} ${SRC}/$(CHECK_ORDER).cpp ${VSIIE_Solver}.o  $(CCLINKLIBS) -o $(CHECK_ORDER) 

${VSIIE_Solver}.o:  ${DIR_VSIIE}/VSIIE_Solver.cpp   
	$(CC) -c  $(C_FLAGS)  -I${DIR_IVP} -I${INC} ${DIR_VSIIE}/VSIIE_Solver.cpp  -o  ${VSIIE_Solver}.o 


.PHONY: clean

clean:
	rm  $(OBJ_SOLVER)  *.o 
cleanout:
	rm  *.txt 
