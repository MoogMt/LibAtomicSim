module conversion

	! Description
	!---------------------------------------------------------------------
	! Contains various conversions variables in order to allow conversion between 
	! units easily without looking them up all the time
	!---------------------------------------------------------------------

	! Geometry
	!---------------------------------------------------------------------
	double precision, parameter :: pi = 3.14159265359
	double precision, parameter :: degr2rad = pi/180
	double precision, parameter :: rad2degr = 180/pi
	!---------------------------------------------------------------------

	! Pressure 
	!---------------------------------------------------------------------
	double precision, parameter :: aU2Pascal = 2.9421912E13
	double precision, parameter :: pascal2AU = 1/aU2Pascal
	double precision, parameter :: gpa2au = 3.399E-5
	double precision, parameter :: au2Gpa = 1/gpa2au
	double precision, parameter :: kbar2Gpa = 0.1
	double precision, parameter :: GPa2kbar = 10
	!---------------------------------------------------------------------

	! Frequency
	!---------------------------------------------------------------------
	double precision, parameter :: tHz2cm = 33.3565
	double precision, parameter :: cm2THz = 1/tHz2cm
	!---------------------------------------------------------------------

	! Energy
	!---------------------------------------------------------------------
	double precision, parameter :: hartree2SI = 4.35974417E18;
	double precision, parameter :: sI2hartree = 1/hartree2SI
	double precision, parameter :: ry2eV = 13.605693009
	double precision, parameter :: eV2Ry = 1/ry2eV
	double precision, parameter :: kelvin2hartree = 0.0000032
	double precision, parameter :: hartree2kelvin = 1/kelvin2hartree
	!---------------------------------------------------------------------

	! Frequency
	!---------------------------------------------------------------------
	double precision, parameter :: hatime2fs = 2.4188843265857*0.01
	double precision, parameter :: fs2hatime = 1/hatime2fs
	double precision, parameter :: fs2ps = 0.001
	double precision, parameter :: ps2fs = 0.001
	!---------------------------------------------------------------------

	! Distance
	!---------------------------------------------------------------------
	double precision, parameter :: bohr2Ang = 0.5291772109
	double precision, parameter :: ang2Bohr = 1/bohr2Ang
	double precision, parameter :: meters2Ang = 10E-10
	double precision, parameter :: meters2nm  = 10E-9
	double precision, parameter :: nm2meters = 1/meters2nm
	double precision, parameter :: ang2meters = 1/meters2Ang
	!---------------------------------------------------------------------

	! Volume
	!---------------------------------------------------------------------
	double precision, parameter :: a3tocm3 = 1E-24
	double precision, parameter :: a3tom3  = 1E-30
	!---------------------------------------------------------------------

	! Mass
	!---------------------------------------------------------------------
	double precision, parameter :: amu2gram = 1.660539040E-24
	double precision, parameter :: gram2Amu = 1/amu2gram
	!---------------------------------------------------------------------

end module conversion