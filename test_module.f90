program test_module
  use xyz
  integer::n_steps=5, n_atoms=96

  double precision,dimension(:,:,:),allocatable :: traj

  integer:: size_len=len_trim("/home/moogmt/Data/CO2/CO2_AIMD/8.82/3000K/TRAJEC.xyz")
  call readxyztraj("/home/moogmt/Data/CO2/CO2_AIMD/8.82/3000K/TRAJEC.xyz",size_len,n_steps,n_atoms,traj)
  write(*,*) "test"

  write(*,*) traj(n_steps,n_atoms,1), traj(n_steps,n_atoms,2), traj(n_steps,n_atoms,3)
end program test_module
