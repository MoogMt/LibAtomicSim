module utils

character(len=240), parameter :: all_names = "H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr" &
												// "RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeLaCePrNdPmSmEuGdTbDyHoErTmYb"   &
												// "LuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"                                   &
												// "FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOg"

contains

! Check that the file existence
!---------------------------------------------------------------------
! Not that useful given that inquire exists already...
subroutine check_existence( file_path, found )
	implicit none

	character*200, intent(in) :: file_path
	logical, intent(out) :: found

	! Check whether file exists at filepath
	inquire( file=file_path, exist=found )

	! If the file is not found, sends a message and returns false
	if ( .not. found ) then
		write(*,*) "file not found", file_path
	! Else returns true
	else
		found=.true.
	endif 

	return
end subroutine check_existence
!---------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine z2Names( z, atom_name )
	! All intent have to be declared
	implicit none

	! Argument
	integer, intent(in) :: z

	! Output 
	character(len=2), intent(out) :: atom_name


	return
end subroutine z2Names
!---------------------------------------------------------------------

end