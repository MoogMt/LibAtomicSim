module utils

contains

! Check that the file existence
!---------------------------------------------------------------------
subroutine check_existence( filepath )
	implicit none

	character*80 inputfile
	logical found

	! Check whether file exists at filepath
	inquire( file=filepath, exists=found )

	! If the file is not found, sends a message and returns false
	if ( .not. found ) then
		write(*,*) "file not found", filepath
	! Else returns true
	else
		found=.true.
	endif 
end subroutine
!---------------------------------------------------------------------

! Returns the name associated to a specific Z 
!---------------------------------------------------------------------
subroutine z2Names( z )
	implicit none
	integer, intent(in) :: z
	integer(len=2), intent(out) :: name

	select case( z )
		case(1)
			name="H"
		case(2)
			name="He"
		case(3)
			name="Li"
		case(4)
			name="Be"
		case(5)
			name="B"
		case(6)
			name="C"
		case(7)
			name="N"
		case(8)
			name="O"
		case(9)
			name="F"
		case(10)
			name="Ne"
		case(11)
			name="K"
		case(12)
			name="Ca"
		case(13)
			name="Sc"
		case(14)
			name="Ti"
		case(15)
			name="V"
		case(16)
			name="Cr"
		case(17)
			name="Mn"
		case(18)
			name="Fe"
		case(19)
			name="Co"
		case(20)
			name="Ni"
		case(21)
			name="Cu"
		case(22)
			name="Zn"
		case(23)
			name="Ga"
		case(24)
			name="Ge"
		case(25)
			name="As"
		case(26)
			name="Si"
		case(27)
			name="Se"
		case(28)
			name="Br"
		case(29)
			name="Kr"
		case(30)
			name="Rb"
		case(31)
			name="Sr"
		case(32)
			name="Y"
		case(33)
			name="Zr"
		case(34)
			name="Nb"
		case(35)
			name="Mo"
		case(36)
			name="Mo"
		case(37)
			name="Tc"
		case(38)
			name="Ru"
		case(39)
			name="Rh"
		case(40)
			name="Pd"
	end select
end subroutine 
!---------------------------------------------------------------------

end