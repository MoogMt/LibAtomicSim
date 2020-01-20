program xyztomolecules
!
implicit none
!
double precision, allocatable :: x(:,:)
integer, allocatable :: a(:,:),typ(:),used(:),try(:),nextry(:)
double precision :: box(3), cutoff2(4,4), dx,dd
character  :: symb*2, line*200, filename*32, arg*32
character*2 :: styp(4)
integer :: i,j,k,h,n,l,step, histotyp(4),ntotmol,nmol(10000),mol(10000,4)
logical :: newmol

!!!!!!!!!!!!!!!!! INPUT
if (iargc().lt.4) then
  write(*,*) "ERROR: usage is:     xyz2molecules.x filename a b c"
  stop
endif
call getarg(1, filename)
call getarg(2, arg)
read(arg,*) box(1)
call getarg(3, arg)
read(arg,*) box(2)
call getarg(4, arg)
read(arg,*) box(3)
!!!!!!!!!!!!!!!!! SETUP
l=0
step=0
cutoff2(1,1)=1.0D0**2  ! H-H
do i=2,4
  cutoff2(1,i)=1.2D0**2  ! heavy-H
  cutoff2(i,1)=1.2D0**2  ! heavy-H
  do j=2,4
    cutoff2(i,j)=1.7D0**2  ! heavy-heavy
  enddo
enddo
styp(1)="H"
styp(2)="C"
styp(3)="O"
styp(4)="N"
!!!!!!!!!!!!!!!!! PARSE
open(30,file=filename,status="old")
do
  read(30,'(a200)',end=333) line
  l=l+1
  if (l.eq.1) then
    i=0
    step=step+1
    if (step.eq.1) then
      read(line,*) n
      write (*,*) "n ",n
      allocate(x(n,3),a(n,n),typ(n),used(n),try(n),nextry(n))
    endif
    write (*,'(A,I6,4X,$)') "STEP ",step
  endif
  if (l.gt.2) then
    i=i+1
!    write (*,*) "i ",i
    read(line,*) symb,x(i,1:3)
    if     (symb.eq.'H') then
      typ(i)=1
    elseif (symb.eq.'C') then
      typ(i)=2
    elseif (symb.eq.'O') then
      typ(i)=3
    elseif (symb.eq.'N') then
      typ(i)=4
    else
      write(*,*) "ERROR: unknown type of atom ",symb
    endif
    if (i.eq.n) then
      l=0
      !!!!!!!!!!!!!!!! ANALYSIS
      ! build adjacency matrix
      a(:,:)=0
      do j=1,n-1
        do k=j+1,n
          dd=0.D0
          do h=1,3
            dx=x(j,h)-x(k,h)
            if (dx> box(h)*0.5) dx=dx-box(h)
            if (dx<-box(h)*0.5) dx=dx+box(h)
            dd=dd+dx*dx
          enddo
          if (dd.le.cutoff2(typ(j),typ(k))) then
            a(j,k)=1
            a(k,j)=1
          endif
        enddo
      enddo
!      do j=1,n
!        write(55,'(1000i2)') a(j,:) ! DEBUG
!      enddo
      ! find molecules
      ntotmol=0
      used(:)=0
      do j=1,n
        if (used(j).eq.1) cycle
        ! start new molecule
        histotyp(:)=0
        try(:)=0
        try(j)=1
        do
          ! start new series of neighbors
          nextry(:)=0
          do k=1,n
            if (try(k).eq.0) cycle
            used(k)=1
            histotyp(typ(k))=histotyp(typ(k))+1
            do h=1,n
              if ((used(h).eq.0).and.(a(k,h).eq.1)) nextry(h)=1
            enddo
          enddo
          try(:)=nextry(:)
          if (sum(try(:)).eq.0) exit
        enddo
        ! store molecule
        newmol=.true.
        do h=1,ntotmol
          if (sum((histotyp(:)-mol(h,:))**2).eq.0) then
            newmol=.false.
            nmol(h)=nmol(h)+1
          endif
        enddo
        if (newmol) then
          ntotmol=ntotmol+1
          nmol(ntotmol)=1
          mol(ntotmol,:)=histotyp(:)
        endif
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! print molecules
      do h=1,ntotmol
        write(*,'(I4,A,$)') nmol(h)," "
        do k=1,4
          if (mol(h,k).gt.0 .and.mol(h,k).lt.10 ) write(*,'(A,I1,$)') trim(styp(k)),mol(h,k)
          if (mol(h,k).ge.10.and.mol(h,k).lt.100) write(*,'(A,I2,$)') trim(styp(k)),mol(h,k)
          if (mol(h,k).ge.100)                    write(*,'(A,I3,$)') trim(styp(k)),mol(h,k)
        enddo
        write(*,'(4X,$)')
      enddo
      write(*,*) ""
    endif
    ! new step
  endif
enddo
!
333 continue
!
end program xyztomolecules
