real(8) function get_theta(x1,x2,y1,y2) result(theta_out)
  real(8), intent(in) :: x1,x2,y1,y2
  real(8) :: topofeq, botofeq
  topofeq = y2-y1
  botofeq = x2-x1
  if (topofeq.eq.0..or.botofeq.eq.0.) then
     theta_out = 0.
  else
     theta_out = atan(topofeq/botofeq)
  endif
end function get_theta

real(8) function get_phi(x1,x2,y1,y2,z1,z2) result(phi_out)
  real(8), intent(in) :: x1,y1,x2,y2,z1,z2
  phi_out = atan(sqrt((x2-x1)**2.+(y2-y1)**2.) / (z2-z1))
end function get_phi

program test
  use SHTOOLS
  use OMP_LIB
  implicit none

  integer :: i,j,k,l,n,o,iatm,jatm,err,m,m1,m2,m3,lmin,lmax,katm, &
       & n_neighbours,mm,nc,jj,tnc,oo,identity,filenumber
  integer :: ntime, natoms, step, size, comp,ncomp,latm,matm,matmm
  integer, dimension(:), allocatable :: clist,checkoff,members
  integer, dimension(:,:), allocatable :: tracked_cluster

  real(8) :: rando,a,b,c,r,dx,pi, theta,phi,rq,rq2,solid_criteria, &
       & bondlength,get_phi,get_theta,rlim,Q6
  real(8), dimension(66) :: p
  real(8), dimension(13) :: w3j
  real(8), dimension(:), allocatable :: box,q4,what6,w6
  real(8), dimension(:,:), allocatable :: xyz,qbar4m,qbar6m,dotq6m,q6tilda,solids
  real(8), dimension(:,:,:), allocatable :: neighbours

  character(len=14) :: char_timestep
  character(len=11) :: filename
  character(len=6) :: deleteline
  character(len=49) :: catcommand,sed_remove
  logical :: busy,found_step,checklist

  !
  ! We Now require the sub script has this: seq $(grep 'ITEM: TIMESTEP' pos.lammps | wc -l)
  ! Perhaps with a statement to avoid if checklist already exists
  !

  allocate(box(3))
  pi = 3.14159265359
  solid_criteria = 0.7
  bondlength = 2.89
  rlim = bondlength
  rq = bondlength*1.1
  rq2 = rq**2.


  busy = .true.
  do while (busy.eqv..true.)
     inquire(file='pos.lammps',opened=busy)
     if (busy.eqv..false.) then
        ntime = 0
        open(1,file='pos.lammps',status='old',action='read')
        i = 0
        do
           read(1,'(A)',end=100,iostat=err) char_timestep
           if (err.ne.0) then
              busy = .true.
              goto 100
           endif
           if (char_timestep.eq.'ITEM: TIMESTEP') ntime = ntime + 1
           i = i + 1
        enddo
100     continue
        rewind(1)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) natoms
        close(1)
     endif
  enddo
  print *, ntime, natoms

  allocate(clist(ntime)) ! 1 = not started, 0 = started, 2 = completed
  allocate(xyz(natoms,4))
  clist = 0


  inquire(file='checklist',exist=checklist)
  if (checklist.eqv..true.) then
     busy = .true.
     do while (busy.eqv..true.)
        inquire(file='checklist',opened=busy)
        if (busy.eqv..false.) then
           open(2,file='checklist',status='old',action='read')
           do
              read(2,*,end=200,iostat=err) i
              if (err.ne.0) then
                 busy = .true.
                 goto 200
              endif
              clist(i) = 1
           enddo
200        continue
           close(2)
        endif
     enddo
  endif

  do while (sum(clist(:)).gt.0)
     found_step = .false.


     inquire(file='checklist',exist=checklist)
     if (checklist.eqv..true.) then
        busy = .true.
        do while (busy.eqv..true.)
           inquire(file='checklist',opened=busy)
           if (busy.eqv..false.) then
              open(2,file='checklist',status='old',action='read',iostat=err)
              if (err.ne.0) then
                 busy = .true.
              else
                 do
                    read(2,*,end=201,iostat=err) i
                    if (err.ne.0) then
                       busy = .true.
                       goto 201
                    endif
                    clist(i) = 1
                 enddo
201              continue
                 close(2)
              endif
           endif
        enddo
     endif





     do while (found_step.eqv..false.)
        call random_number(rando)
        i = int(rando*real(ntime+1))
        if (i.ge.1.and.i.le.ntime) then
           step = i
           if (clist(step).eq.1) then
              found_step = .true.
              clist(step) = 0
           endif
        endif
     enddo
     write(sed_remove,fmt='(a9,i0,a13)') "sed -i '/",step,"/d' checklist"
     call execute_command_line(sed_remove)

     !print *,step
     xyz = 0.

     found_step = .false.
     j = 0
300  continue
     do while (found_step.eqv..false.)
        inquire(file='pos.lammps',opened=busy)
        if (busy.eqv..false.) then
           open(1,file='pos.lammps',status='old',action='read')
           do i = 1, step-1
              do j = 1, natoms+9
                 read(1,*)
              enddo
           enddo
           do i = 1, 9
              read(1,*)
           enddo
           do i = 1, natoms
              read(1,*,iostat=err) j,k,a,b,c
              if (err.ne.0) then
                 busy = .true.
                 close(1)
                 goto 300
              endif
              xyz(j,1) = a
              xyz(j,2) = b
              xyz(j,3) = c
              xyz(j,4) = real(k)
           enddo
           close(1)
           found_step = .true.
        endif
     enddo

     ! meat and potatoes
     print *, 'Getting into it now'
     allocate(q4(natoms))
     allocate(w6(natoms))
     allocate(what6(natoms))

     allocate(qbar4m(natoms,11))
     allocate(qbar6m(natoms,15))
     allocate(dotq6m(natoms,30))
     allocate(q6tilda(natoms,15))

     allocate(tracked_cluster(natoms,natoms+3))
     allocate(checkoff(natoms))
     allocate(solids(natoms,4))
     allocate(members(natoms))
     allocate(neighbours(natoms,30,2))

     q4 = 0.
     what6 = 0.
     qbar4m = 0.
     qbar6m = 0.
     w6 = 0.
     Q6 = 0.
     q6tilda = 0.
     solids = 0.
     members = 0


     do iatm = 1, natoms ! primary atom loop
        n_neighbours = 0
        do jatm = 1, natoms ! neighbour atom loop
           if (iatm.eq.jatm) goto 301
           r = 0.
           do k = 1, 3
              dx = xyz(jatm,k) - xyz(iatm,k)
              if (abs(dx).gt.rq*3.) goto 301
              if (abs(dx).gt.0.5*box(k)) dx = (-dx)**0. * box(k) + dx
              r = r + dx**2.
           enddo
           if (r.gt.rq2) goto 301
           r = sqrt(r)
           do j = 1, 30
              if (neighbours(iatm,j,1).eq.0) then
                 neighbours(iatm,j,1) = real(jatm)
                 neighbours(iatm,j,2) = r
                 n_neighbours = n_neighbours+1
                 goto 302
              endif
           enddo
           a = 0.
           do j = 1, 30
              if (neighbours(iatm,j,2).gt.a) then
                 a = neighbours(iatm,j,2)
                 k = j
              endif
           end do
           if (r.lt.a) then
              neighbours(iatm,k,2) = a
              neighbours(iatm,k,1) = jatm
           endif
302        continue

           theta = get_theta(xyz(iatm,1),xyz(jatm,1),xyz(iatm,2),xyz(jatm,2))
           phi = get_phi(xyz(iatm,1),xyz(jatm,1),xyz(iatm,2),xyz(jatm,2),xyz(iatm,3),xyz(jatm,3))
           call PlmBar(p,4,cos(theta))
           do m = -4, 4
              if (m.ge.0) then
                 qbar4m(iatm,m+5) = qbar4m(iatm,m+5) + (p(int(PlmIndex(4,abs(m)))) * cos(real(abs(m))*phi)) * (r-rq)**2.
              else
                 qbar4m(iatm,m+5) = qbar4m(iatm,m+5) + (p(int(PlmIndex(4,abs(m)))) * sin(real(abs(m))*phi)) * (r-rq)**2.
              endif
           enddo
           qbar4m(iatm,10) = qbar4m(iatm,10) + (r-rq)**2.
           qbar4m(iatm,11) = qbar4m(iatm,11) + 1.
           !######### q6 portion     
           call PlmON(p,6,cos(theta))
           do m = -6, 6
              if (m.ge.0) then
                 qbar6m(iatm,m+7) = qbar6m(iatm,m+7) + (p(int(PlmIndex(6,abs(m)))) * cos(real(abs(m))*phi)) * (r-rq)**2.
              else
                 qbar6m(iatm,m+7) = qbar6m(iatm,m+7) + (p(int(PlmIndex(6,abs(m)))) * cos(real(abs(m))*phi)) * (r-rq)**2.
              endif
           enddo
           qbar6m(iatm,14) = qbar6m(iatm,14) + (r-rq)**2.
           qbar6m(iatm,15) = qbar6m(iatm,15) + 1.
301        continue
        enddo !~ Neighbour atom loop
        do m = -4, 4
           qbar4m(iatm,m+5) = qbar4m(iatm,m+5) / real(n_neighbours)
        enddo
        do m = -6, 6
           qbar6m(iatm,m+7) = qbar6m(iatm,m+7) / real(n_neighbours)
        enddo
        !###### q4 portion
        do m = -4, 4
           q4(iatm) = q4(iatm) + qbar4m(iatm,m+5)**2.
        enddo
        q4(iatm) = sqrt(q4(iatm) * (( 4.*pi)/(2.*4.+1.)))
        !###### what6 portion
        do m = -6, 6
           what6(iatm) = what6(iatm) + qbar6m(iatm,m+7)**2.
        enddo
        !###### w6 portion
        w6(iatm) = 0.
        lmin = 6
        lmax = 6
        do m1 = -6, 6
           do m2 = -6, 6
              do m3 = -6, 6
                 if (m1+m2+m3.eq.0) then                
                    call Wigner3j(w3j,lmin,lmax,6,6,m1,m2,m3)
                    w6(iatm) = w6(iatm) + (w3j(6-abs(m1)+1)* qbar6m(iatm,m1+7)* qbar6m(iatm,m2+7)* qbar6m(iatm,m3+7))
                 endif
              enddo
           enddo
        enddo
        what6(iatm) = w6(iatm) / abs(what6(iatm))**(3./2.)
        if (what6(iatm).ne.what6(iatm)) what6(iatm) = 0.
        if (q4(iatm).ne.q4(iatm)) q4(iatm) = 0.
     enddo!~ Primary atom loop
     !#### Q6 portion
     do m = -6, 6
        Q6 = Q6 + (sum(qbar6m(:,m+7)) / sum(qbar6m(:,14)))**2.
     enddo
     Q6 = sqrt(((4.*pi)/13.) * Q6)
     do iatm = 1, natoms
        a = 0.
        do m = -6, 6
           a = a + qbar6m(iatm,m+7)**2.      
        enddo
        a = sqrt(a)
        do m = -6, 6
           q6tilda(iatm,m+7) = qbar6m(iatm,m+7) / a
        enddo
     enddo
     dotq6m = 0.   
     do iatm = 1, natoms !# Solid ID loop
        k = 0
        l = 0
        a = 0.
        do j = 1, 30
           if (neighbours(iatm,j,1).ne.0.) then
              k = k + 1
              do m = -6, 6
                 dotq6m(iatm,j) = dotq6m(iatm,j) + q6tilda(iatm,m+7) * q6tilda(int(neighbours(iatm,j,1)),m+7)
              enddo
              a = a + dotq6m(iatm,j)
              if (dotq6m(iatm,j).ge.solid_criteria.and.neighbours(iatm,j,2).le.rlim) then
                 l = l + 1
              endif
           endif
        enddo

        if (l.ge.8) then
           solids(iatm,1) = 1.
           solids(iatm,2) = what6(iatm)
           solids(iatm,3) = q4(iatm)
        endif
     enddo  !# Solid ID loop
     do iatm = 1, natoms
        k = 0
        do j = 1, 30
           if (neighbours(iatm,j,1).ne.0..and.neighbours(iatm,j,2).lt.bondlength) then
              jatm = int(neighbours(iatm,j,1))
              if (solids(iatm,1).eq.1.and.solids(jatm,1).eq.1) k = k + 1
           endif
        enddo
        if (k.ge.8)  members(iatm) = 1
     enddo
     checkoff = 0
     mm = 0
     nc = 0
     tracked_cluster = 0

     do iatm = 1, natoms
        if (members(iatm).eq.1) then
           if (checkoff(iatm).eq.0) then
              checkoff(iatm) = 1
              m = 1
              nc = nc + 1
              tnc = tnc + 1
              tracked_cluster(iatm,m) = iatm
              n = m
1311          continue
              l = n ! go from last max
              n = m ! rest to current max
              do o = l, n
                 jatm = tracked_cluster(iatm,o)
                 do j = 1, 30
                    if (neighbours(jatm,j,1).ne.0..and.neighbours(jatm,j,2).lt.bondlength) then
                       katm = int(neighbours(jatm,j,1))
                       if (checkoff(katm).eq.0.and.members(katm).eq.1) then
                          m = m + 1
                          tracked_cluster(iatm,m) = katm
                          checkoff(katm) = 1
                       endif
                    endif
                 enddo
              enddo
              if (m.gt.n) goto 1311            
              identity = identity + 1
              do o = 1, m
                 do oo = 1, m
                    katm = tracked_cluster(iatm,o)
                    jatm = tracked_cluster(iatm,oo)
                    if (katm.eq.jatm) then
                       theta = get_theta(xyz(katm,1),xyz(jatm,1),xyz(katm,2),xyz(jatm,2))
                       phi = get_phi(xyz(katm,1),xyz(jatm,1),xyz(katm,2),xyz(jatm,2),xyz(katm,3),xyz(jatm,3))
                    endif
                 enddo
              enddo
              jj = jj + count(tracked_cluster(iatm,:).gt.0)

              ! Compositional part
              comp = 0
              ncomp = 0
              do latm = 1, natoms
                 do matm = 1, count(tracked_cluster(iatm,:).ne.0)
                    ! stop checking if the subject atom is in the cluster
                    if (latm.eq.tracked_cluster(iatm,matm)) goto 5670
                 enddo
                 !now check the distances between this subject atom and all present in the cluster
                 do matmm = 1, count(tracked_cluster(iatm,:).ne.0)
                    matm = tracked_cluster(iatm,matmm)
                    if (matm.ne.0.) then
                       r = 0.
                       do k = 1, 3
                          dx = xyz(matm,k) - xyz(latm,k)
                          if (abs(dx).gt.bondlength*2.) goto 5671 ! if this atom doesnt work, move to next in cluster
                          if (abs(dx).gt.0.5*box(k)) dx = (-dx)**0. * box(k) + dx
                          r = r + dx**2.
                       enddo
                       r = r**0.5
                       if (r.lt.bondlength) then
                          ncomp = ncomp + 1
                          if (xyz(latm,4).eq.1) comp = comp + 1
                          goto 5670 ! if any atoms in the cluster are within bonding distance, this atom qualifies so we can stop checking
                       endif
                    endif
5671                continue
                 enddo
5670             continue
              enddo
              tracked_cluster(iatm,natoms+2) = comp
              tracked_cluster(iatm,natoms+3) = ncomp
           endif
        endif
     enddo
     o = count(solids(:,4).ne.0.)
     do iatm = 1, natoms
        katm = count(tracked_cluster(iatm,:).ne.0)
        if (katm.gt.0) then
           m = 0
           do jatm = 1, katm
              if (xyz(tracked_cluster(iatm,jatm),4).eq.1) m = m + 1
           enddo

           filenumber = int(rando*100000)
           write(filename,fmt='(a5,i6.6)') 'fort.',filenumber
           open(2040,file=filename,status='replace')
           write(2040,*) step, katm, real(m)/real(katm), tracked_cluster(iatm,natoms+2), tracked_cluster(iatm,natoms+3) !,tracked_cluster(iatm,:)
           close(2040)
           write(catcommand,fmt='(a4,a11,a23,a11)') 'cat ',filename,' >> clusters.out && rm ',filename    
           call execute_command_line(catcommand)
           INQUIRE(FILE='checklist', SIZE=size)
           if (size.eq.0) then
              print *, "EXITING"
              call exit(0)
           endif
        endif
     enddo
     clist(step) = 0
     deallocate(q4)
     deallocate(w6)
     deallocate(what6)
     deallocate(qbar4m)
     deallocate(qbar6m)
     deallocate(dotq6m)
     deallocate(q6tilda)
     deallocate(tracked_cluster)
     deallocate(checkoff)
     deallocate(solids)
     deallocate(members)
     deallocate(neighbours)

  enddo
  !print *, "DONE"



end program test
