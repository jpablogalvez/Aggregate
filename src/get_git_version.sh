#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)

cat <<EOF > version.f90
!======================================================================!
!
      module version
!
      implicit none
!
      contains
!
!======================================================================!
!
      subroutine print_version
!
      write(*,'(2X,A)') "Commit id    :  $git_hash"
      write(*,'(2X,A)') "Commit date  :  $git_date"
!          
      return
      end subroutine print_version
!
!======================================================================!
!
      end module version
!
!======================================================================!
EOF
