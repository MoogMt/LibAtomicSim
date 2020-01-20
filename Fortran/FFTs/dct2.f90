!-------
! DCT2
!-------------------------------------------------------------------------------------------------
! Performs the discrete cosine transform of type II on velocity autocorrelation to get the vdos
!-------------------------------------------------------------------------------------------------
!http://en.wikipedia.org/wiki/Discrete_cosine_transform DCT-II
program dct2
  
  implicit none
  
  !functions
  integer::get_N_lines,argcount
  
  ! From Hz to cm-1
  real*8,parameter::thz2cm1=33.3565d0
  
  integer:: n,tmp_int
  real*8 :: pif, frame, norm,total_time,tmp_dble
  integer :: i,j,k,nat,nindex,kk
  integer,allocatable,dimension(:)::atom_index
  !real*8 :: x(0:n-1),y(0:n-1),xnew(0:n-1),z(0:n-1),spec_all(0:n-1),t(0:n-1)
  real*8,allocatable,dimension(:) :: t
  real*8,allocatable,dimension(:) :: x_all,y_all,xnew_all,z_all,modulus_all
  real*8,allocatable,dimension(:) :: x_spec,y_spec,xnew_spec,z_spec,modulus_spec
  real*8,allocatable,dimension(:,:) :: x,y,xnew,z,modulus
  character(len=200)::wq_char,tmp_char,index_file,vv_file='vv.dat',prefix='spectrum'
  logical::spec_on=.false.
  
  ! Getting flags and arguments
  argcount=IARGC()
  
  ! If no flags or arguments, prints help
  if(argcount==0) then
     ! Usage
     print'(a)', 'USAGE: /!\ SINCE WE WENT FROM DCT TO DCT2 THE FREQ ARE PROBABLY OFF'
     print'(a)', './dct2.x -N_atoms n_atoms -index_file file.index -vv_file vv.dat -out prefix'
     print'(a)', '[-N_atoms n_atoms       ] NECESSARY - the numer of atoms '
     print'(a)', '[-index_file file.index ] OPTIONNAL - index (start at 0) of atoms whose contribution is to be computed as the sum of the individual VACF - format: N_index\n i1 i2 ... iN on the same line'
     print'(a)', '[-vv_file vv.dat        ] OPTIONNAL - file where the VACF is stored (format: t a1(t) a2(t) ... aN(t) - t in fs - default=vv.dat'
     print'(a)', '[-out prefix            ] OPTIONNAL - prefix for the spectrum, program adds .real'
     stop
  endif

  ! GETTING FLAGS ARGUMENTS
  do i=1,argcount
     call GETARG(i,wq_char)
     ! Necessary, get the number of atoms
     IF(INDEX(wq_char,'-N_atoms').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)nat
        CYCLE
     ! Optionnal, name of input file
     ELSE IF(INDEX(wq_char,'-vv_file').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)vv_file
        CYCLE
     ! Optionnal, name of output file
     ELSE IF(INDEX(wq_char,'-out').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)prefix
        CYCLE
     ! Optionnal, index file to compute only some atoms
     ELSE IF(INDEX(wq_char,'-index_file').NE.0)THEN
        spec_on=.true.
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)index_file
        CYCLE
     ENDIF
  enddo
  

  ! If index file, reads the index to count as one
  if(spec_on) then
     open(50,file=index_file,status='old')
     read(50,*),nindex
     allocate(atom_index(nindex))
     read(50,*)(atom_index(i),i=1,nindex)
     close(50)
  else
     nindex=0
  endif
  
  ! read data
  open(40,file=vv_file,status='old')
  ! Get number of line of file
  n=get_N_lines(40)
  ! File pointer back to the begining of file
  rewind(40)
  ! Allocate tables with expected number of points
  allocate(t(n))  ! Time (in fs)
  allocate(x(nat,n),x_all(n),x_spec(n)) 
  allocate(y(nat,n),y_all(n),y_spec(n)) 
  allocate(xnew(nat,n),xnew_all(n),xnew_spec(n))
  allocate(z(nat,n),z_all(n),z_spec(n))
  allocate(modulus(nat,n),modulus_all(n),modulus_spec(n))

  ! Reads the velocity autocorrelation
  do i=1,n
     ! Reading line i
     read(40,*) t(i),(x(k,i),k=1,nat)

     ! Summing all contributions from all atoms
     tmp_dble=0.d0
     do k=1,nat
        tmp_dble = tmp_dble + x(k,i)
     enddo
     ! Final correlation
     x_all(i)=tmp_dble

     ! Summing contributions from specified atom indexes
     tmp_dble=0.d0
     do k=1,nindex
        kk=atom_index(k)+1 ! +1 because they are index (they start at 0)
        tmp_dble=tmp_dble+x(kk,i)
     enddo
     ! Final correlation
     x_spec(i)=tmp_dble
     
  end do
  ! Close the file
  close(40)

  ! Computes the total time of simulation (assumes unit fs)
  total_time=(t(n)-t(1))*1.d-3

  ! Computation constant
  pif=4*datan(1.d0)/n

  ! Initialize vector
  y=0.d0
  z=0.d0
  y_all=0.d0
  z_all=0.d0
  y_spec=0.d0
  z_spec=0.d0

  ! DCT Type 2
  ! - Loop to get the components
  do i=1,n
     ! Loop to get the sum
     do j=1,n 
        ! Compute the FFT component
        tmp_dble=dcos(pif*(j-0.5)*(i-1))
        do k=1,nat
           y(k,i)=y(k,i)+x(k,j)*tmp_dble
        enddo
        ! FFT for all atoms
        y_all(i)=y_all(i)+x_all(j)*tmp_dble
        ! FFT for only specific atoms
        y_spec(i)=y_spec(i)+x_spec(j)*tmp_dble
     enddo
  enddo

  ! Normalization
  !y=y/dble(n)
  !y_all=y_all/dble(n)
  !y_spec=y_spec/dble(n)
  
  ! Open output file for all atoms
  open(unit=32, file=trim(adjustl(prefix))//'.real', action = 'write')
  ! Write DCT2 to file for all atoms
  do i=1,n
     write(32,'(2f20.12,$)') DBLE(i-1)/2.d0*1./total_time*thz2cm1, y_all(i)
     do k=1,nat
        write(32,'(f20.12,$)') y(k,i)
     enddo
     write(32,'(a)') ' '
  end do
  close(32)
  ! Output file for specified atoms
  if(spec_on) then
     open(unit=62, file=trim(adjustl(prefix))//'_index.real', action = 'write')
     do i=1,n
        write(62,'(3f20.12)') DBLE(i-1)/2.d0*1./total_time*thz2cm1, y_all(i),y_spec(i)
     end do
     close(62)
  endif
  
  ! backward discrete cosine transform (the inverse of DCT-II is DCT-III*2/n)
  !pif=4*datan(1.d0)/n
  
  !do k=1,nat
  !   xnew(k,:)=0.5d0*y(k,1)
  !enddo
  !xnew_all(:)=0.5d0*y_all(1)
  !xnew_spec(:)=0.5d0*y_spec(1)
  
 ! do i=1,n
  !   do j=2,n
  !      tmp_dble=dcos(pif*(j-1.d0)*(i-0.5d0))
  !      do k=1,nat
  !         xnew(k,i)=xnew(k,i)+y(k,j)*tmp_dble
  !      enddo
  !      xnew_all(i)=xnew_all(i)+y_all(j)*tmp_dble
  !      xnew_spec(i)=xnew_spec(i)+y_spec(j)*tmp_dble
  !   enddo
  !enddo
  !xnew=xnew*2.d0/n
  !xnew_all=xnew_all*2.d0/n
  !xnew_spec=xnew_spec*2.d0/n
  ! 
  !do i=1,n
  !   do k=1,nat
  !      if(abs(x(k,i)-xnew(k,i))>=0.0000001d0) write(*,*) x(k,i), xnew(k,i), x(k,i)-xnew(k,i)
  !   enddo
  !   if(abs(x_all(i)-xnew_all(i))>=0.0000001d0) write(*,*) x_all(i), xnew_all(i), x_all(i)-xnew_all(i)
  !   if(abs(x_spec(i)-xnew_spec(i))>=0.0000001d0) write(*,*) x_spec(i), xnew_spec(i), x_spec(i)-xnew_spec(i)
  !end do
  
end program dct2



!----------------------------------------------------------------------!
! Get the number of lines in a file
!----------------------------------------------------------------------!
function get_N_lines(unit_number)
  implicit none
  
  integer::unit_number,get_N_lines
  integer::ios = 0, N_lines = 0
  
  rewind(unit_number)
  
  do while(ios == 0)
     read(unit_number,*,iostat=ios)
     if(ios == 0) N_lines = N_lines +1
  enddo
  
  if(N_lines == 0) then
     write(*,'(a,i4)') 'ERROR : non-existing or empty file in unit :', unit_number
     stop
  endif
  
  get_N_lines = N_lines
  
  rewind(unit_number)
end function get_N_lines
!----------------------------------------------------------------------!
