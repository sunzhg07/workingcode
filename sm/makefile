LIBS = -lblas -llapack -lm
FV = gfortran  -g  -fbounds-check -Wall
all = sm_module.o sm_basis.o

sm.exe :	${all}
	   ${FV} -o sm.exe ${all} ${LIBS}

%.o: %.f90
	${FV} -c $< ${LIBS}

clean:
	rm *.mod *.exe *.o *__genmod.f90
		
