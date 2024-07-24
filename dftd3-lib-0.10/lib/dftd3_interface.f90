! Interface for c& c++
! AR. Allouche
subroutine dftd3c(instr, version, nAtoms, atnum, ccoords, edisp, cgrads)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*), intent(IN) :: instr
  integer, parameter :: wp = kind(1.0d0)
  integer :: i,j,k
  integer :: nAtoms
  integer :: version
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)

  real(wp) :: stress(3,3)
  real(wp), allocatable :: grads(:,:)
  real(wp), allocatable :: coords(:,:)
  character(len=100) :: func=" "
  integer :: l

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  l=0
  do
      if (instr(l+1) == C_NULL_CHAR) exit
      func(l+1:l+1) = instr(l+1)
      l = l + 1
  end do

  allocate(grads(3,nAtoms))
  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  !call dftd3_set_functional(dftd3, func='dftb3', version=version, tz=.false.)
  call dftd3_set_functional(dftd3, func=func, version=version, tz=.false.)

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, edisp, grads)
  
  k=1
  do i=1,3
  do j=1,nAtoms
     cgrads(k) = grads(i,j)
     k = k + 1
  enddo
  enddo
  deallocate(grads)
  deallocate(coords)

end subroutine dftd3c

subroutine dftd3cpbc(instr, version, nAtoms, atnum, ccoords, clatVecs,  edisp, cgrads, cstress)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*), intent(IN) :: instr
  integer, parameter :: wp = kind(1.0d0)
  integer :: i,j,k
  integer :: nAtoms
  integer :: version
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)
  real(wp) :: clatVecs(9)
  real(wp) :: cstress(9)

  real(wp) :: stress(3,3)
  real(wp) :: latVecs(3,3)
  real(wp), allocatable :: grads(:,:)
  real(wp), allocatable :: coords(:,:)
  character(len=100) :: func=" "
  integer :: l

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  l=0
  do
      if (instr(l+1) == C_NULL_CHAR) exit
      func(l+1:l+1) = instr(l+1)
      l = l + 1
  end do

  allocate(grads(3,nAtoms))
  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  k=1
  do i=1,3
  do j=1,3
     latVecs(i,j) = clatVecs(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  call dftd3_set_functional(dftd3, func=func, version=version, tz=.false.)

  ! Calculate dispersion and gradients for periodic case
  call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, edisp, grads, stress)
  k=1
  do i=1,3
  do j=1,3
     cstress(k) = stress(i,j)
     k = k + 1
  enddo
  enddo
  k=1
  do i=1,3
  do j=1,nAtoms
     cgrads(k) = grads(i,j)
     k = k + 1
  enddo
  enddo
  deallocate(grads)
  deallocate(coords)

end subroutine dftd3cpbc

! Without gradients
subroutine dftd3ce(instr, version, nAtoms, atnum, ccoords, edisp)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*), intent(IN) :: instr
  integer, parameter :: wp = kind(1.0d0)
  integer :: i,j,k
  integer :: nAtoms
  integer :: version
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)

  real(wp), allocatable :: coords(:,:)
  character(len=100) :: func=" "
  integer :: l

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  l=0
  do
      if (instr(l+1) == C_NULL_CHAR) exit
      func(l+1:l+1) = instr(l+1)
      l = l + 1
  end do

  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  !call dftd3_set_functional(dftd3, func='dftb3', version=version, tz=.false.)
  call dftd3_set_functional(dftd3, func=func, version=version, tz=.false.)

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, edisp)
  
  deallocate(coords)

end subroutine dftd3ce

subroutine dftd3cpbce(instr, version, nAtoms, atnum, ccoords, clatVecs,  edisp)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*), intent(IN) :: instr
  integer, parameter :: wp = kind(1.0d0)
  integer :: i,j,k
  integer :: nAtoms
  integer :: version
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)
  real(wp) :: clatVecs(9)

  real(wp) :: latVecs(3,3)
  real(wp), allocatable :: coords(:,:)
  character(len=100) :: func=" "
  integer :: l

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  l=0
  do
      if (instr(l+1) == C_NULL_CHAR) exit
      func(l+1:l+1) = instr(l+1)
      l = l + 1
  end do

  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  k=1
  do i=1,3
  do j=1,3
     latVecs(i,j) = clatVecs(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  call dftd3_set_functional(dftd3, func=func, version=version, tz=.false.)

  ! Calculate dispersion and gradients for periodic case
  call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, edisp)
  deallocate(coords)

end subroutine dftd3cpbce

! Interface for c& c++
! AR. Allouche
! pars : # Parameters are as follows: (s6, rs6, s18, rs18, alp6)
!        # otherwise known as         (s6,  a1,  s8, a2, alp6)

subroutine dftd3cpars(pars, version, nAtoms, atnum, ccoords, edisp, cgrads)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  real(wp) :: pars(5)
  integer, intent(in) :: version
  integer :: i,j,k
  integer :: nAtoms
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)

  real(wp) :: stress(3,3)
  real(wp), allocatable :: grads(:,:)
  real(wp), allocatable :: coords(:,:)

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  allocate(grads(3,nAtoms))
  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  call dftd3_set_params(dftd3, pars=pars, version=version)

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, edisp, grads)
  
  k=1
  do i=1,3
  do j=1,nAtoms
     cgrads(k) = grads(i,j)
     k = k + 1
  enddo
  enddo
  deallocate(grads)
  deallocate(coords)

end subroutine dftd3cpars

subroutine dftd3cpbcpars(pars, version, nAtoms, atnum, ccoords, clatVecs,  edisp, cgrads, cstress)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  real(wp) :: pars(5)
  integer, intent(in) :: version
  integer :: i,j,k
  integer :: nAtoms
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)
  real(wp) :: clatVecs(9)
  real(wp) :: cstress(9)

  real(wp) :: stress(3,3)
  real(wp) :: latVecs(3,3)
  real(wp), allocatable :: grads(:,:)
  real(wp), allocatable :: coords(:,:)

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  allocate(grads(3,nAtoms))
  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo


  k=1
  do i=1,3
  do j=1,3
     latVecs(i,j) = clatVecs(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  call dftd3_set_params(dftd3, pars=pars, version=version)

  ! Calculate dispersion and gradients for periodic case
  call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, edisp, grads, stress)
  k=1
  do i=1,3
  do j=1,3
     cstress(k) = stress(i,j)
     k = k + 1
  enddo
  enddo
  k=1
  do i=1,3
  do j=1,nAtoms
     cgrads(k) = grads(i,j)
     k = k + 1
  enddo
  enddo
  deallocate(grads)
  deallocate(coords)

end subroutine dftd3cpbcpars

! Without gradients
subroutine dftd3cepars(pars, version, nAtoms, atnum, ccoords, edisp)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  real(wp) :: pars(5)
  integer, intent(in) :: version
  integer :: i,j,k
  integer :: nAtoms
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)

  real(wp), allocatable :: coords(:,:)

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  call dftd3_set_params(dftd3, pars=pars, version=version)

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, edisp)
  
  deallocate(coords)

end subroutine dftd3cepars

subroutine dftd3cpbcepars(pars, version, nAtoms, atnum, ccoords, clatVecs,  edisp)
  use dftd3_api
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  real(wp) :: pars(5)
  integer, intent(in) :: version
  integer :: i,j,k
  integer :: nAtoms
  integer :: atnum(nAtoms)
  real(wp) :: ccoords(3*nAtoms)
  real(wp) :: cgrads(3*nAtoms)
  real(wp) :: clatVecs(9)

  real(wp) :: latVecs(3,3)
  real(wp), allocatable :: coords(:,:)

  ! Same conversion factor as in dftd3
  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  real(wp) :: edisp

  allocate(coords(3,nAtoms))
  k=1
  do i=1,3
  do j=1,nAtoms
     coords(i,j) = ccoords(k)
     k = k + 1
  enddo
  enddo

  k=1
  do i=1,3
  do j=1,3
     latVecs(i,j) = clatVecs(k)
     k = k + 1
  enddo
  enddo

  ! Initialize dftd3
  call dftd3_init(dftd3, input)


  call dftd3_set_params(dftd3, pars=pars, version=version)

  ! Calculate dispersion and gradients for periodic case
  call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, edisp)
  deallocate(coords)

end subroutine dftd3cpbcepars
