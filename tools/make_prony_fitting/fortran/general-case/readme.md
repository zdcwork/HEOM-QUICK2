you should first change the variable 'NAME' in Makefile to the correct path of the current directory,
and then use command [make] to generate the executable file 'prony.x'.
run [./prony.x > out &] and we have
eta    -->   cb.data
eta^dag  -->   cd.data
gamma  -->   cgama.data

chang and save environmental parameters in the ./main/main.f90 files and command [make] the directory again.

the program in '../general-case' works for arbitrary spectral function.
the program in '../Lorentzian' works for the Lorentzian-type spectral function.
