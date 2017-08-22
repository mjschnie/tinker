c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readhorton  --  get information from Horton     ##
c     ##                             ISA output                      ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readhorton" reads in the Horton ISA output in cartesian global 
c     multipoles
c
c
      subroutine readhorton
      use sizes
      use atomid
      use atoms
      use dma
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k
      integer idma,next
      integer freeunit
      real*8 term
      logical exist,done
      logical use_bohr
      character*3 atmnam
      character*120 record
      character*120 dmafile
      real*8 globalpole(20)
      character*5 shape
      integer length
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(mp))  allocate (mp(maxatm))
      if (.not. allocated(dpx))  allocate (dpx(maxatm))
      if (.not. allocated(dpy))  allocate (dpy(maxatm))
      if (.not. allocated(dpz))  allocate (dpz(maxatm))
      if (.not. allocated(q20))  allocate (q20(maxatm))
      if (.not. allocated(q21c))  allocate (q21c(maxatm))
      if (.not. allocated(q21s))  allocate (q21s(maxatm))
      if (.not. allocated(q22c))  allocate (q22c(maxatm))
      if (.not. allocated(q22s))  allocate (q22s(maxatm))
c
c     zero out the atomic coordinates and DMA values
c
      print *,"made it in"
      do i = 1, maxatm
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         mp(i) = 0.0d0
         dpx(i) = 0.0d0
         dpy(i) = 0.0d0
         dpz(i) = 0.0d0
         q20(i) = 0.0d0
         q21c(i) = 0.0d0
         q21s(i) = 0.0d0
         q22c(i) = 0.0d0
         q22s(i) = 0.0d0
      end do
c
c     try to get a filename from the command line arguments
c
      call nextarg (dmafile,exist)
      if (exist) then
         call basefile (dmafile)
         call suffix (dmafile,'horton','old')
         inquire (file=dmafile,exist=exist)
      end if
c
c     ask for the user specified Horton output filename
c
      do while (.not. exist)
         write (iout,10)
 10      format (/,' Enter Horton ISA Output File Name :  ',$)
         read (input,20)  dmafile
 20      format (a120)
         call basefile (dmafile)
         call suffix (dmafile,'dma','old')
         inquire (file=dmafile,exist=exist)
      end do
      print *,"got file"
c
c     first open and then read the Horton output file
c
      idma = freeunit ()
      open (unit=idma,file=dmafile,status='old')
c
c     get coordinates and multipoles from Horton output file
c     HOW DO I GET COORDINATES???
c
      i = 0
      rewind (unit=idma)
      do while (.true.)
         read (idma,30,err=50,end=50)  record
 30      format (a120)
         if (record(1:28) .eq. 'Dataset,cartesian_multipoles') then
c
c     find number of atoms
c
            read (idma,*,err=50,end=50) shape,n,length
c
c     perform dynamic allocation of some global arrays
c
            if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
c
c     read in global multipoles
c
            do i = 1, n
               read (idma,*,err=50,end=50) globalpole(:)
               print *,"globalpole",globalpole
               do k = 1, 13
                  rpole(k,i) = globalpole(k)
               end do
            end do
         else if (record(1:14) .eq. 'Dataset,change') then
            goto 50
         end if
      end do
 50   continue
      n = i
c
c     convert coordinates from Bohrs to Angstroms if needed
c
c     not sure what units horton is using...
c
c     find atomic numbers in verbose GDMA output if available
c
c     
      done = .false.
      rewind (unit=idma)
      do while (.true.)
         read (idma,80,err=100,end=100)  record
 80          format (a120)
         if (record(1:16) .eq. 'Nuclear charges:') then
            k = min(n,20)
            read (record(17:120),*,err=100,end=100)  (atomic(i),i=1,k)
            do while (k .ne. n)
               j = k + 1
               k = min(n,k+20)
               read (idma,90,err=100,end=100)  record
 90                      format (a120)
               read (record,*,err=100,end=100)  (atomic(i),i=j,k)
            end do
            done = .true.
         end if
      end do
 100   continue
      close (unit=idma)
c
c     attempt to get atomic numbers from GDMA atom names
c
      if (.not. done) then
         do i = 1, n
            atomic(i) = 0
            atmnam = name(i)
            call upcase (atmnam)
            if (atmnam(1:2) .eq. 'SI') then
               atomic(i) = 14
            else if (atmnam(1:2) .eq. 'CL') then
               atomic(i) = 17
            else if (atmnam(1:2) .eq. 'BR') then
               atomic(i) = 35
            else if (atmnam(1:1) .eq. 'H') then
               atomic(i) = 1
            else if (atmnam(1:1) .eq. 'B') then
               atomic(i) = 5
            else if (atmnam(1:1) .eq. 'C') then
               atomic(i) = 6
            else if (atmnam(1:1) .eq. 'N') then
               atomic(i) = 7
            else if (atmnam(1:1) .eq. 'O') then
               atomic(i) = 8
            else if (atmnam(1:1) .eq. 'F') then
               atomic(i) = 9
            else if (atmnam(1:1) .eq. 'P') then
               atomic(i) = 15
            else if (atmnam(1:1) .eq. 'S') then
               atomic(i) = 16
            else if (atmnam(1:1) .eq. 'I') then
               atomic(i) = 53
            else
               read (atmnam,*,err=110,end=110)  atomic(i)
 110                     continue
            end if
         end do
      end if
c
c     print the global frame Cartesian atomic multipoles
c
      write (iout,120)
 120  format (/,' Global Frame Cartesian Multipole Moments :')
      do i = 1, n
         write (iout,130)  i,name(i),atomic(i)
 130     format (/,' Site:',i8,9x,'Name:',3x,a3,7x,'Atomic Number:',i8)
         write (iout,140)  x(i),y(i),z(i)
 140     format (/,' Coordinates:',5x,3f15.6)
         write (iout,150)  rpole(1,i)
 150     format (/,' Charge:',10x,f15.5)
         write (iout,160)  rpole(2,i),rpole(3,i),rpole(4,i)
 160     format (' Dipole:',10x,3f15.5)
         write (iout,170)  rpole(5,i)
 170     format (' Quadrupole:',6x,f15.5)
         write (iout,180)  rpole(8,i),rpole(9,i)
 180     format (18x,2f15.5)
         write (iout,190)  rpole(11,i),rpole(12,i),rpole(13,i)
 190     format (18x,3f15.5)
      end do
c     
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for use as traceless values
c
c     keep this if horton gives moments in bohr like gdma
      do i = 1, n
         do k = 2, 4
            rpole(k,i) = rpole(k,i) * bohr
         end do
         do k = 5, 13
            rpole(k,i) = rpole(k,i) * bohr**2 / 3.0d0
         end do
      end do
      return
      end
c
c
c
