program su2mshscale

implicit none

integer             :: i,j,k, ierr,ndim,ielem
integer             :: elem_type,elem_line(15)
integer             :: Nelem, Npoin, point_start,point_end,line_count=0,elem_start
real*8              :: x1,y1,z1,scaling=0.001
real*8              :: x2,y2,z2
character(len=200)  :: fname,text,word,fout
character(len=20)   :: ndim_str,nelem_str,npoin_str

fname = 'AltTEtopology_blade_full.su2'
fout = 'AltTEtopology_blade_scaled_full.su2'
ndim_str = 'NDIME='
nelem_str = 'NELEM='
npoin_str = 'NPOIN='
open(unit=10, file=fname(1:len_trim(fname)),status='old')
open(unit=11, file=fout(1:len_trim(fout)),status='unknown')

do

read (10,"(a)",iostat=ierr) text ! read line into character variable
line_count = line_count + 1
if (ierr /= 0) exit
read (text,*) word ! read first word of line

!First look for nDim
   if (word == ndim_str) then ! found search string at beginning of line
      read (text,*) word,ndim
      print*,"ndim =",ndim
   end if
! Now look for Nelem
   if (word == nelem_str) then ! found search string at beginning of line
      read (text,*) word,Nelem
      elem_start=line_count+1
      print*,"nelem =",Nelem,elem_start
   end if
! Now look for Npoin and find out which line point starts from
   if (word == npoin_str) then ! found search string at beginning of line
      read (text,*) word,Npoin
      point_start = line_count+1
      point_end = point_start + Npoin - 1
      print*,"npoin =",Npoin,point_start,point_end
   end if

enddo

! allocate(x2(Npoin),y2(Npoin),z2(Npoin))

rewind(10)
! Write dimension
write(11,'(a6,i3)') ndim_str(1:len_trim(ndim_str)),ndim

line_count = 0
do i=1,elem_start-1
read(10,*)
enddo

! Write number of elements
write(11,'(a6,i9)') nelem_str(1:len_trim(nelem_str)),nelem
write(*,'(i9)')nelem
!Read and write elements into new file
do i=elem_start,elem_start+Nelem-1
read(10,"(a)",iostat=ierr) text
! print*,text
! read(text,*) elem_line

write(11,'(a)') text(1:len_trim(text))
enddo

!Commented lines
do i=elem_start+Nelem,point_start-1
read(10,*)
enddo

! Start reading points
write(11,'(a6,i9)') npoin_str,Npoin
write(*,'(i9)')Npoin
do i=point_start,point_end
if (nDim .eq. 2) then
read(10,*)x1,y1,ielem
x2 = x1*scaling
y2 = y1*scaling
write(11,*) x2,y2,ielem
else
read(10,*)x1,y1,z1,ielem
x2 = x1*scaling
y2 = y1*scaling
z2 = z1*scaling
write(11,*) x2,y2,z2,ielem
endif
enddo

do
read (10,"(a)",iostat=ierr) text ! read line into character variable
if (ierr /= 0) exit
write(11,'(a)')text(1:len_trim(text))
enddo

if (nDim .eq. 2) then
print*,x1,y1,ielem
else
print*,x1,y1,z1,ielem
endif




close(10)
close(11)

! deallocate(x2,y2,z2)
end program su2mshscale
