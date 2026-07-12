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
       real(kind=8),dimension(:),allocatable           ::  xagg     !  Aggregate-count distribution
       real(kind=8),dimension(:),allocatable           ::  pagg     !  Molecule-weighted aggregate distribution
       real(kind=8),dimension(:,:),allocatable         ::  pqagg    !  Species-conditioned aggregate distribution
       real(kind=8),dimension(:),allocatable           ::  xsize    !  Aggregate-count size distribution
       real(kind=8),dimension(:),allocatable           ::  psize    !  Molecule-weighted size distribution
       real(kind=8),dimension(:,:),allocatable         ::  pqsize   !  Species-conditioned size distribution
       real(kind=8),dimension(:,:),allocatable         ::  homxsize !  Homogeneous aggregate-count size distribution by species
       real(kind=8),dimension(:),allocatable           ::  mixxsize !  Mixed aggregate-count size distribution
       real(kind=8),dimension(:,:),allocatable         ::  hompsize !  Homogeneous molecule-weighted size distribution by species
       real(kind=8),dimension(:),allocatable           ::  mixpsize !  Mixed molecule-weighted size distribution
       real(kind=8),dimension(:,:),allocatable         ::  hompop   !  Homogeneous aggregate populations by species and size
       real(kind=8),dimension(:,:),allocatable         ::  homconc  !  Homogeneous aggregate concentrations by species and size
       real(kind=8),dimension(:,:),allocatable         ::  homnum   !  Homogeneous aggregate counts by species and size
       real(kind=8),dimension(:),allocatable           ::  mixpop   !  Mixed aggregate populations by size
       real(kind=8),dimension(:),allocatable           ::  mixconc  !  Mixed aggregate concentrations by size
       real(kind=8),dimension(:),allocatable           ::  mixnum   !  Mixed aggregate counts by size
       real(kind=8),dimension(:),allocatable           ::  mixapop  !  Mixed aggregate populations by identifier
       real(kind=8),dimension(:),allocatable           ::  mixaconc !  Mixed aggregate concentrations by identifier
       real(kind=8),dimension(:),allocatable           ::  mixanum  !  Mixed aggregate counts by identifier
       real(kind=8)                                    ::  sumagg   !  Total number of aggregate samples
       real(kind=8)                                    ::  summol   !  Total number of molecular samples
       real(kind=8),dimension(:),allocatable           ::  sumqmol  !  Total molecular samples by species
       integer                                         ::  msize    !  Maximum aggregate size
       integer                                         ::  nmax     !  Maximum aggregate identifier
!
       end module properties
!
!======================================================================!
