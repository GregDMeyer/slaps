
CFLAGS 	        =
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
EXAMPLESC       = petsc_impl.c

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

petsc_impl: petsc_impl.o  chkopts
	-${CLINKER}  -o petsc_impl petsc_impl.o ${PETSC_MAT_LIB}
	${RM} petsc_impl.o

include ${PETSC_DIR}/lib/petsc/conf/test
