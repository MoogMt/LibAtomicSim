program prout
	
	use omp_lib

INTEGER :: thread_id

!$OMP PARALLEL PRIVATE(thread_id)
    thread_id = OMP_GET_THREAD_NUM()
    DO i=0,OMP_GET_MAX_THREADS()
        IF (i == thread_id) THEN
            PRINT *, "Hello from process: ", thread_id
        END IF
        !$OMP BARRIER
    END DO
!$OMP END PARALLEL

end program prout

