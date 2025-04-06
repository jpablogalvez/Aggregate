!======================================================================!
!
       module screening
!
       use omp_var,  only:  np,chunkscrn
!
       use omp_lib
!
       implicit none
!
       contains
!
!======================================================================!
!
! SCRNINT - SCReeNing INTeractions
!
! This subroutine determines if any interaction between two molecules is
!  present in the previous and current snapshots or in the current and
!  next snapshots (remove collisions) as well as in the previous and the
!  next snapshots (add oscillations)
!
       subroutine scrnint(nnode,oldadj,adj,newadj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
! Local variables
! 
       integer                                       ::  i,j     !
!
! Comparing the diagonal blocks of two adjacency matrices
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,oldadj,newadj,nnode)                      &
!$omp             private(i,j)                                         &
!$omp             schedule(dynamic,chunkscrn)
!
       do i = 1, nnode-1
         do j = i+1, nnode
           if ( (oldadj(j,i).and.adj(j,i)) .or.                        &
                (adj(j,i).and.newadj(j,i)) .or.                        &
                                    (oldadj(j,i).and.newadj(j,i)) ) then
             adj(j,i) = .TRUE.
             adj(i,j) = .TRUE.
           else
             adj(j,i) = .FALSE.
             adj(i,j) = .FALSE.
           end if
         end do
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine scrnint
!
!======================================================================!
!
! SCRNINT - SCReeNing COLlision
!
! This subroutine determines if any interaction between two molecules is
!  present in the previous and current snapshots or in the current and
!  next snapshots (remove collisions)
!
       subroutine scrncol(nnode,oldadj,adj,newadj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
! Local variables
! 
       integer                                       ::  i,j     !
!
! Comparing the diagonal blocks of two adjacency matrices
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,oldadj,newadj,nnode)                      &
!$omp             private(i,j)                                         &
!$omp             schedule(dynamic,chunkscrn)
!
       do i = 1, nnode-1
         do j = i+1, nnode
           if ( (.NOT.(oldadj(j,i).and.adj(j,i))) .and.                &
                                (.NOT.(adj(j,i).and.newadj(j,i))) ) then
             adj(j,i) = .FALSE.
             adj(i,j) = .FALSE.
           end if
         end do
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine scrncol
!
!======================================================================!
!
! SCRNOSC - SCReeNing OSCillations
!
! This subroutine determines if any interaction between two molecules is
!  present in the previous and the next snapshots (add oscillations)
!
       subroutine scrnosc(nnode,oldadj,adj,newadj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
! Local variables
! 
       integer                                       ::  i,j     !
!
! Comparing the diagonal blocks of two adjacency matrices
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,oldadj,newadj,nnode)                      &
!$omp             private(i,j)                                         &
!$omp             schedule(dynamic,chunkscrn)
!
       do i = 1, nnode-1
         do j = i+1, nnode
           if ( oldadj(j,i) .and. newadj(j,i) ) then
             adj(j,i) = .TRUE.
             adj(i,j) = .TRUE.
           end if
         end do
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine scrnosc
!
!======================================================================!
!
       end module screening
!
!======================================================================!
