module utils

contains

character(len=240) :: all_species_names= "H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr" //&
										"RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeLaCePrNdPmSmEuGdTbDyHoErTmYb" // &
										"LuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn" // &
										"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOg"

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

end