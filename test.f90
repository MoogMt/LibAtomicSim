PROGRAM TEST

 use MPI

 integer :: nb_procs, rang, code

 call MPI_INIT(code)

 call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
 call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

 write(*,*) "Je suis le processus", rang, "parmi", nb_procs

 call MPI_FINALIZE(code)

END PROGRAM TEST 
