module xyz

contains

! Descriptions
!  - Functions that handle *.xyz files


! Get number of atoms in the traj (assumes a real trajectory, not a set of structures)
!-----------------------------------------------------------------------------------------
subroutine getnbatoms( file_path , str_len, handle_nb, nb_atoms )
  ! All intents have to be explicitely written
  implicit none

  ! Arguments
  integer, intent(in) :: str_len                  ! Length of the string containing the file path
  integer, intent(in) :: handle_nb                ! number of the file handle
  character(len=str_len), intent(in) :: file_path ! String containing the input file name
  ! Output
  integer, intent(out) :: nb_atoms                ! number of atoms 

  ! Opens input file
  open( handle_nb, file=file_path, status="old" )

  ! Get number of atoms as the first element of the first line of the file
  read( handle_nb, * ) nb_atoms

  ! Close file
  close( handle_nb )
end subroutine getnbatoms
!-----------------------------------------------------------------------------------------

! Get number of atoms and steps in the traj (assume a tral trajectory, not to use in a set of structures of different sizes)
!--------------------------------------------------------------
subroutine getnbatomssteps( file_path, str_len, handle_nb, nb_atoms, nb_steps )
  ! All intents have to be explicitely written
  implicit none

  ! Arguments 
  integer, intent(in) :: str_len                  ! length of the string of the path/name of the input file
  integer, intent(in) :: handle_nb                ! handler number for the file handling
  character(len=str_len), intent(in) :: file_path ! name or path of the input file
  ! Output 
  integer, intent(out) :: nb_steps                ! number of steps in the trajectory
  integer, intent(out) :: nb_atoms                ! number of atoms in the trajectory
  ! Local variables
  integer :: sel                                  ! frame counter, counts the lines within a single trajectory frame
  integer :: err                                  ! error handler for the file manipulation

  ! Get number of atoms
  call getnbatoms( file_path, str_len, handle_nb, nb_atoms )

  ! Init output variable
  nb_steps=0

  ! Init frame counter
  sel=1 

  ! Opens input file
  open( handle_nb, file=file_path, status="old" )

  ! Loop over all lines of the file
  do while ( err == 0 )
    ! If the frame counter is at 1, increments step counter 
    if ( sel == 1  ) then
      nb_steps = nb_steps + 1
    ! When at the end of the frame reinitiate frame counter to 0
    elseif ( sel == nb_atoms + 2 ) then
       sel = 0
    endif

    ! Read line, sends status into err
    read( handle_nb, *, iostat=err )

    ! Increments counter at all lines
    sel=sel+1
  enddo

  ! Closes input file
  close( handle_nb )

  ! Step counter is always one over the correct number
  nb_steps = nb_steps - 1
end subroutine getnbatomssteps
!--------------------------------------------------------------

! Read an xyz traj and returns positions
!--------------------------------------------------------------
subroutine readxyztraj( file_path, str_len, handle_nb, n_steps, n_atoms, positions )
  ! All variables intent have to be declared
  implicit none

  ! Arguments
  integer, intent(in) :: str_len
  integer, intent(in) :: handle_nb
  character(len=str_len), intent(in):: file_path

  ! Output
  integer, intent(out) :: n_steps
  integer, intent(out) :: n_atoms
  double precision, dimension(:,:,:),allocatable,intent(out):: positions

  ! Local variables
  character(len=4) :: dummy_name
  integer:: atom, step

  ! Get number of atoms and steps in the trajectory
  call getnbatomssteps( file_path, str_len, handle_nb, n_atoms, n_steps )

  ! Allocate tensor for atom positions in traj
  allocate( positions(3, n_atoms, n_steps) )

  ! Opens input file
  open( handle_nb, file=file_path, status="old" )
  do step=1,n_steps
    ! Read the number of atom line
    read( handle_nb, * )

    ! Read the comment line
    read( handle_nb, * )

    ! Loop over atoms
    do atom=1,n_atoms
      ! Read line for each atom with atom name and atom positions 
      read( handle_nb, * ) dummy_name, positions( 1, atom, step ), positions( 2, atom, step ), positions( 3, atom, step )
    enddo
  enddo

  ! Close input file
  close( handle_nb )
end subroutine readxyztraj
!--------------------------------------------------------------

end module xyz
