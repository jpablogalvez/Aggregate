!======================================================================!
!
       module properties
!
       implicit none
!
! Average properties
!
       real(kind=8),dimension(:,:,:),allocatable       ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:),allocatable           ::  pop      !  Populations
       real(kind=8),dimension(:),allocatable           ::  conc     !  Concentrations
       real(kind=8),dimension(:),allocatable           ::  prob     !  Concentrations
       real(kind=8),dimension(:),allocatable           ::  frac     !  Molar fractions
       real(kind=8),dimension(:),allocatable           ::  cin      !  Stechiometric concentration
       real(kind=8)                                    ::  volu     !  Simulation box volume
       integer,dimension(:),allocatable                ::  num      !  Number of molecules
       integer                                         ::  msize    !  Maximum aggregate size
       integer                                         ::  nmax     !  Maximum aggregate identifier
!
       end module properties
!
!======================================================================!
