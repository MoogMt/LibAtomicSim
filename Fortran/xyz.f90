module xyz

use utils

contains

! Description of module
!---------------------------------------------------------------------
! Contains functions that handle *.xyz files
!---------------------------------------------------------------------

! Get number of atoms in the traj (assumes a real trajectory, not a set of structures)
!-----------------------------------------------------------------------------------------
subroutine getnbatoms( file_path , str_len, handle_nb, nb_atoms )
  ! All intents have to be explicitely written
  implicit none

  !----------------------------------------
  ! Arguments
  integer, intent(in) :: str_len                  ! Length of the string containing the file path
  integer, intent(in) :: handle_nb                ! number of the file handle
  character(len=str_len), intent(in) :: file_path ! String containing the input file name
  ! Output
  integer, intent(out) :: nb_atoms                ! number of atoms 
  !----------------------------------------

  ! Opens input file
  open( handle_nb, file=file_path, status="old" )

  ! Get number of atoms as the first element of the first line of the file
  read( handle_nb, * ) nb_atoms

  ! Close file
  close( handle_nb )

  ! Ends subroutine
  return
end subroutine getnbatoms
!-----------------------------------------------------------------------------------------
! Get number of atoms and steps in the traj (assume a tral trajectory, not to use in a set of structures of different sizes)
!--------------------------------------------------------------
subroutine getnbatomssteps( file_path, str_len, handle_nb, nb_atoms, nb_steps )
  ! All intents have to be explicitely written
  implicit none

  !----------------------------------------
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
  integer :: number_line                          ! number of lines in the file
  !----------------------------------------

  ! Get number of atoms
  call getnbatoms( file_path, str_len, handle_nb, nb_atoms )

  ! Get number of lines in the file
  call getNumberLine( file_path, str_len, handle_nb, number_line )

  ! Computes number of steps
  nb_steps = number_line/(nb_atoms+2)

  ! Ends subroutine
  return
end subroutine getnbatomssteps
!--------------------------------------------------------------

! Read an xyz traj and returns positions
!--------------------------------------------------------------
subroutine readxyztraj( file_path, str_len, handle_nb, n_steps, n_atoms, positions, atoms_z )
  ! All variables intent have to be declared
  implicit none

  !----------------------------------------
  ! Arguments
  integer, intent(in) :: str_len                 ! Length of the string for path of input file
  integer, intent(in) :: handle_nb               ! handle number
  character(len=str_len), intent(in):: file_path ! path to input file
  ! Output
  integer, intent(out) :: n_steps                                          ! number of steps
  integer, intent(out) :: n_atoms                                          ! number of atoms
  double precision, dimension(:,:,:), allocatable, intent(out):: positions ! positions of atoms
  integer, dimension(:), allocatable, intent(out) :: atoms_z               ! atomic z
  ! Local variables
  character(len=4) :: name ! Dummy variable for atom names
  integer :: z_prime       ! Dummy for atomic number
  integer :: atom          ! Dummy variable for atom loop
  integer :: step          ! Dummy variable for step loop
  !----------------------------------------

  ! Get number of atoms and steps in the trajectory
  call getnbatomssteps( file_path, str_len, handle_nb, n_atoms, n_steps )

  ! Allocate tensor for atom positions in traj
  allocate( positions(3, n_atoms, n_steps) )

  ! Allocate vector for atomic numbers
  allocate( atoms_z(n_atoms) )

  ! Opens input file
  open( handle_nb, file=file_path, status="old" )

  ! Loop over steps
  do step=1,n_steps
    ! Read the number of atom line
    read( handle_nb, * )

    ! Read the comment line
    read( handle_nb, * )

    ! Loop over atoms
    do atom=1,n_atoms
      ! Read line for each atom with atom name and atom positions 
      read( handle_nb, * ) name, positions( 1:3, atom, step )

      ! Use first step to get the atomic numbers
      if ( step .eq. 1 ) then
        call names2Z( name, z_prime )
        atoms_z(atom) = z_prime
      endif

    enddo
  enddo

  ! Close input file
  close( handle_nb )

  ! Ends subroutine
  return
end subroutine readxyztraj
!--------------------------------------------------------------

end module xyz
