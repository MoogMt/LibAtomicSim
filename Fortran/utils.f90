module utils


! Description of module
!---------------------------------------------------------------------
! Contains various functions that are very useful and that did not fit
! elsewhere
!---------------------------------------------------------------------


! General variables
!---------------------------------------------------------------------
! number of elements
integer, parameter :: max_z = 120  
! Super String containing all species names
character(len=max_z*2), parameter :: all_names = "H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr" &
												// "RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeLaCePrNdPmSmEuGdTbDyHoErTmYb"       &
												// "LuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"                                       &
												// "FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOg"
!---------------------------------------------------------------------

contains

! File dealing functions
!---------------------------------------------------------------------
! Not that useful given that inquire exists already...
subroutine checkExistence( file_path, file_path_length, found )
	implicit none

	!----------------------------------------
	! Argument
	integer, intent(in) :: file_path_length                  ! length of the string for the path of input file
	character(len=file_path_length), intent(in) :: file_path ! path to the input file 
	! Output
	logical, intent(out) :: found  ! Whether the file was found or not
	!----------------------------------------

	! Check whether file exists at filepath
	inquire( file=file_path, exist=found )

	! If the file is not found, sends a message and returns false
	if ( .not. found ) then
		write(*,*) "file not found", file_path
		found=.false.
	! Else returns true
	else
		found=.true.
	endif 

	! Ends subroutine
	return
end subroutine checkExistence
! Get number of line in a given file
subroutine getNumberLine( file_path, file_length, handle_nb, number_line )
	! All intent have to be explicit
	implicit none

	!----------------------------------------
	! Argument
	integer, intent(in) :: handle_nb                    ! Number of the handle of the file
	integer, intent(in) :: file_length                  ! Size of the string of the path of input file
	character(len=file_length), intent(in) :: file_path ! Path to the input file
	! Output
	integer, intent(out) :: number_line ! Number of line in the file
	! Local Variable
	integer :: error_status             ! Reading status of the file
	!----------------------------------------

	! Initialize number of line counter
	number_line=0

	! Opens input file
	open( handle_nb, file=file_path, status="old" )

	! Loop over all lines of file
	do while ( error_status == 0 )
		! Reading line 
		read( handle_nb, *, iostat=error_status )

		! Increments counter
		number_line = number_line + 1
	enddo

	! Decrementing because of first step
	number_line = number_line - 1 

	! Close input file
	close( handle_nb )

	! End subroutine
	return
end subroutine
!---------------------------------------------------------------------

! Transform name and atomic number of atoms
!---------------------------------------------------------------------
! Returns name of an atom from its atomic number
subroutine z2Names( z, atom_name )
	! All intent have to be declared
	implicit none

	!----------------------------------------
	! Argument
	integer, intent(in) :: z                   ! Atomic number
	! Output 
	character(len=2), intent(out) :: atom_name ! name of the target atom
	!----------------------------------------

	! Check that the atomic number if neither too high or too low
	if ( ( z .lt. max_z ) .and. ( z .gt. 0 ) ) then
		! Return name of atoms
		atom_name = all_names( 2*(z-1)+1 : 2*(z-1)+2 )
	else
		! If there is a problem, put "XX" as the atom name
		atom_name = "XX"
	endif

	! Ends subroutine
	return
end subroutine z2Names
! Returns atomic number of an atom from its name
subroutine names2Z( atom_name, z )
	! All intent have to be declared
	implicit none

	!----------------------------------------
	! Argument
	character(len=2), intent(in) :: atom_name ! Name of the atom
	! Output 
	integer, intent(out) :: z                 ! Atomic number
	! Local argument
	integer :: z_prime                        ! Dummy variable for loop
	!----------------------------------------

	! Loop over atomic numbers
	do z_prime=1,max_z
		! IF the atom name matches the one of the databse, returns 
		! the Z of the loop
		if ( atom_name == all_names( 2*(z_prime-1)+1 : 2*(z_prime-1)+2 ) ) then
			z = z_prime
			return
		endif
	enddo

	! If we reach this point we send an error message
	write(*,*) "Atom name was not found in database"

	! Put 0 as the name of the atom
	z=0

	return
end subroutine names2Z
!---------------------------------------------------------------------

end