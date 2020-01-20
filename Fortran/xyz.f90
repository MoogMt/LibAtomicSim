module xyz

contains
! ------------------------------------------------------------------------------
! Read and write XYZ files
! ------------------------------------------------------------------------------

subroutine getnbatoms( file_path , str_len, nb_atoms )

  implicit none

  integer, intent(in) :: str_len
  character(len=str_len), intent(in) :: file_path

  integer, intent(out) :: nb_atoms

  open(25,file=file_path,status="old")
  read(25,*) nb_atoms
  close(25)

end subroutine getnbatoms

subroutine getnbatomssteps( file_path, str_len, nb_atoms, nb_steps )

  implicit none

  integer, intent(in) :: str_len
  character(len=str_len), intent(in) :: file_path

  integer :: err

  integer, intent(out) :: nb_steps
  integer, intent(out) :: nb_atoms
  integer :: sel

  sel=1
  nb_atoms=1
  nb_steps=0

  open(25,file=file_path,status="old",iostat=err)
  read(25,*) nb_atoms ! Getting nb of atoms
  rewind(25)
  do while ( err == 0 )
    if ( sel == 1  ) then
      nb_steps = nb_steps + 1
    elseif ( sel == nb_atoms + 2 ) then
       sel = 0
    endif
    read(25,*,iostat=err)
    sel=sel+1
  enddo
  close(25)

  nb_steps=nb_steps-1

  write(*,*) nb_atoms, nb_steps

end subroutine getnbatomssteps

subroutine readxyztraj( file_path, str_len, n_steps, n_atoms, traj )

  implicit none

  integer,intent(in):: str_len
  character(len=str_len), intent(in):: file_path

  integer:: atom, step

  integer, intent(out) :: n_steps, n_atoms
  character(len=4) :: dummy_name

  double precision, dimension(:,:,:),allocatable,intent(out):: traj

  call getnbatomssteps(file_path,str_len,n_atoms,n_steps)

  allocate( traj(n_steps, n_atoms, 3) )

  open(24,file=file_path,status="old")
  do step=1,n_steps
    ! Read the two dummy lines
    read(24,*)
    read(24,*)
    ! Read one step worth of trajectory
    do atom=1,n_atoms
      read(24,*) dummy_name, traj(step,atom,1),traj(step,atom,2),traj(step,atom,3)
    enddo
  enddo
  close(24)

end subroutine readxyztraj

end module xyz
