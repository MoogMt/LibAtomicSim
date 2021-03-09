program prout
	
	use utils
	use xyz
	use omp_lib

	integer :: z
	character(len=2):: atom_name
	INTEGER :: thread_id

	character(len=200)::path_file
	integer :: n_atoms
	integer :: n_steps

	integer :: handle_xyz
	integer :: str_len

	integer :: number_line

	integer ::i 

	double precision, dimension(:,:,:), allocatable :: positions
	integer, dimension(:), allocatable :: atoms_z

	str_len= 200

	handle_xyz=20

	path_file="/Volumes/B2O3_QUENCH/Project_Quench/NPT/2000-1998/traj.xyz"

	call readtrajxyz( path_file, str_len, handle_xyz, n_steps, n_atoms, positions, atoms_z )

	do i=1,n_atoms
		write(*,*) i, atoms_z(i), positions(:,i,1)
	enddo	


! !$OMP PARALLEL PRIVATE(thread_id)
!     thread_id = OMP_GET_THREAD_NUM()
!     DO i=0,OMP_GET_MAX_THREADS()
!         IF (i == thread_id) THEN
!             PRINT *, "Hello from process: ", thread_id
!         END IF
!         !$OMP BARRIER
!     END DO
! !$OMP END PARALLEL

end program prout

