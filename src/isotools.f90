!======================================================================!
!
       module isotools
!
       use lengths,  only:  lenout
       use printings, only: print_end
!
       implicit none
!
       private
!
       integer,parameter  ::  leniso = 32
       integer,parameter  ::  initial_iso_class_capacity = 64
       integer,parameter  ::  initial_iso_assignment_capacity = 1024
       integer            ::  max_iso_class_capacity = 1000000
!
       type iso_graph
!        Labeled graph used for aggregate-conformation classification.
!        The adjacency matrix is boolean because aggregate representations
!        do not store edge multiplicities, but node labels are chemically
!        meaningful group/body names and are preserved by isomorphism.
         integer                                      ::  n = 0
         integer                                      ::  n_edges = 0
         logical                                      ::  directed = .FALSE.
         logical,dimension(:,:),allocatable           ::  adj
         character(len=leniso),dimension(:),allocatable ::  labels
         integer,dimension(:),allocatable             ::  degree
         integer,dimension(:),allocatable             ::  indegree
         integer,dimension(:),allocatable             ::  outdegree
       end type iso_graph
!
       type iso_class
         type(iso_graph)                              ::  representative
         integer                                      ::  count = 0
         integer                                      ::  old_id = 0
       end type iso_class
!
       type iso_bucket
         logical                                      ::  active = .FALSE.
         logical                                      ::  directed = .FALSE.
         integer                                      ::  iagg = 0
         integer                                      ::  nobs = 0
         integer                                      ::  nclasses = 0
         integer                                      ::  class_capacity = 0
         integer                                      ::  assignment_capacity = 0
         character(len=256)                           ::  labels = ''
         character(len=lenout)                        ::  filebase = ''
         type(iso_class),dimension(:),allocatable     ::  classes
         integer,dimension(:),allocatable             ::  assignments
       end type iso_bucket
!
       type(iso_bucket),dimension(:),allocatable      ::  buckets
!
       public  ::  init_iso,clear_iso,classify_iso_graph,write_iso_results
!
       contains
!
!======================================================================!
!
       subroutine init_iso(nmax,max_iso_classes)
!
       integer,intent(in)  ::  nmax
       integer,intent(in)  ::  max_iso_classes
       integer             ::  i
!
       if ( allocated(buckets) ) call clear_iso()
       max_iso_class_capacity = max_iso_classes
       allocate(buckets(nmax))
       do i = 1, nmax
         buckets(i)%iagg = i
         buckets(i)%active = .FALSE.
         buckets(i)%nobs = 0
         buckets(i)%nclasses = 0
         buckets(i)%class_capacity = 0
         buckets(i)%assignment_capacity = 0
       end do
!
       return
       end subroutine init_iso
!
!======================================================================!
!
       subroutine clear_iso()
!
       integer  ::  i
!
       if ( .not. allocated(buckets) ) return
       do i = 1, size(buckets)
         if ( allocated(buckets(i)%classes) ) deallocate(buckets(i)%classes)
         if ( allocated(buckets(i)%assignments) ) deallocate(buckets(i)%assignments)
       end do
       deallocate(buckets)
!
       return
       end subroutine clear_iso
!
!======================================================================!
!
       subroutine classify_iso_graph(iagg,labels,filebase,adj,directed)
!
       integer,intent(in)                         ::  iagg
       character(len=*),intent(in)                ::  labels
       character(len=*),intent(in)                ::  filebase
       logical,dimension(:,:),intent(in)          ::  adj
       logical,intent(in)                         ::  directed
!
       type(iso_graph)                            ::  graph
       integer                                    ::  class_id
       character(len=lenout)                     ::  base
       integer                                    ::  i
!
       if ( .not. allocated(buckets) ) return
       if ( iagg .lt. 1 .or. iagg .gt. size(buckets) ) return
!
       if ( .not. buckets(iagg)%active ) then
         buckets(iagg)%active = .TRUE.
         buckets(iagg)%directed = directed
         buckets(iagg)%labels = labels
         call iso_filebase(filebase,base)
         buckets(iagg)%filebase = base
         buckets(iagg)%class_capacity = min(initial_iso_class_capacity,&
                                            max_iso_class_capacity)
         buckets(iagg)%assignment_capacity =                           &
              initial_iso_assignment_capacity
         allocate(buckets(iagg)%classes(buckets(iagg)%class_capacity))
         allocate(buckets(iagg)%assignments(                            &
              buckets(iagg)%assignment_capacity))
       end if
!
       call make_graph(labels,adj,directed,graph)
!
       class_id = 0
       do i = 1, buckets(iagg)%nclasses
         if ( are_isomorphic(graph,buckets(iagg)%classes(i)%representative) ) then
           buckets(iagg)%classes(i)%count = buckets(iagg)%classes(i)%count + 1
           class_id = i
           exit
         end if
       end do
!
       if ( class_id .eq. 0 ) then
         call add_iso_class(buckets(iagg),graph,class_id)
       end if
!
       buckets(iagg)%nobs = buckets(iagg)%nobs + 1
       call append_assignment(buckets(iagg),class_id)
!
       return
       end subroutine classify_iso_graph
!
!======================================================================!
!
       subroutine write_iso_results()
!
       integer  ::  i
!
       if ( .not. allocated(buckets) ) return
!
       do i = 1, size(buckets)
         if ( buckets(i)%active ) call write_iso_bucket(buckets(i))
       end do
!
       return
       end subroutine write_iso_results
!
!======================================================================!
!
       subroutine make_graph(labels_line,adj,directed,graph)
!
       character(len=*),intent(in)                ::  labels_line
       logical,dimension(:,:),intent(in)          ::  adj
       logical,intent(in)                         ::  directed
       type(iso_graph),intent(out)                ::  graph
!
       integer                                    ::  n
       integer                                    ::  i,j
!
       n = size(adj,1)
       graph%n = n
       graph%directed = directed
!
       allocate(graph%adj(n,n))
       allocate(graph%labels(n))
       allocate(graph%degree(n),graph%indegree(n),graph%outdegree(n))
!
       call parse_labels(labels_line,n,graph%labels)
!
       graph%adj(:,:) = .FALSE.
       if ( directed ) then
         do i = 1, n
           do j = 1, n
             if ( i .ne. j ) graph%adj(i,j) = adj(i,j)
           end do
         end do
       else
         do i = 1, n-1
           do j = i+1, n
             if ( adj(i,j) .or. adj(j,i) ) then
               graph%adj(i,j) = .TRUE.
               graph%adj(j,i) = .TRUE.
             end if
           end do
         end do
       end if
!
       call compute_graph_properties(graph)
!
       return
       end subroutine make_graph
!
!======================================================================!
!
       subroutine parse_labels(line,n,labels)
!
       character(len=*),intent(in)                  ::  line
       integer,intent(in)                           ::  n
       character(len=leniso),dimension(n),intent(out) ::  labels
!
       integer                                      ::  i
       integer                                      ::  ntok
       integer                                      ::  p0,p1
       integer                                      ::  nline
       character(len=32)                            ::  aux
!
       labels(:) = ''
       nline = len_trim(line)
       i = 1
       ntok = 0
       do while ( i .le. nline .and. ntok .lt. n )
         do while ( i .le. nline )
           if ( line(i:i) .ne. ' ' .and. line(i:i) .ne. char(9) ) exit
           i = i + 1
         end do
         if ( i .gt. nline ) exit
         p0 = i
         do while ( i .le. nline )
           if ( line(i:i) .eq. ' ' .or. line(i:i) .eq. char(9) ) exit
           i = i + 1
         end do
         p1 = i - 1
         ntok = ntok + 1
         labels(ntok) = line(p0:min(p1,p0+leniso-1))
       end do
!
       do i = 1, n
         if ( len_trim(labels(i)) .eq. 0 ) then
           write(aux,'("v",I0)') i
           labels(i) = aux
         end if
       end do
!
       return
       end subroutine parse_labels
!
!======================================================================!
!
       subroutine compute_graph_properties(graph)
!
       type(iso_graph),intent(inout)  ::  graph
       integer                        ::  i,j
!
       graph%degree(:) = 0
       graph%indegree(:) = 0
       graph%outdegree(:) = 0
       graph%n_edges = 0
!
       if ( graph%directed ) then
         do i = 1, graph%n
           do j = 1, graph%n
             if ( graph%adj(i,j) ) then
               graph%outdegree(i) = graph%outdegree(i) + 1
               graph%indegree(j) = graph%indegree(j) + 1
               graph%n_edges = graph%n_edges + 1
             end if
           end do
         end do
         graph%degree(:) = graph%indegree(:) + graph%outdegree(:)
       else
         do i = 1, graph%n
           do j = 1, graph%n
             if ( graph%adj(i,j) ) graph%degree(i) = graph%degree(i) + 1
           end do
         end do
         do i = 1, graph%n-1
           do j = i+1, graph%n
             if ( graph%adj(i,j) ) graph%n_edges = graph%n_edges + 1
           end do
         end do
         graph%indegree(:) = graph%degree(:)
         graph%outdegree(:) = graph%degree(:)
       end if
!
       return
       end subroutine compute_graph_properties
!
!======================================================================!
!
       logical function are_isomorphic(a,b) result(iso)
!
       type(iso_graph),intent(in)  ::  a,b
!
       iso = .FALSE.
       if ( a%n .ne. b%n ) return
       if ( a%directed .neqv. b%directed ) return
       if ( a%n_edges .ne. b%n_edges ) return
       if ( .not. same_tag_multiset(a,b) ) return
!
       if ( a%directed ) then
         iso = are_isomorphic_dir(a,b)
       else
         iso = are_isomorphic_undir(a,b)
       end if
!
       return
       end function are_isomorphic
!
!======================================================================!
!
       logical function are_isomorphic_undir(a,b) result(iso)
!
       type(iso_graph),intent(in)  ::  a,b
       integer                     ::  n
       integer,dimension(:),allocatable  ::  order,map,invmap
       integer,dimension(:),allocatable  ::  candidates
       logical,dimension(:),allocatable  ::  t1,t2,t1tilde,t2tilde
       logical                     ::  found
!
       n = a%n
       allocate(order(n),map(n),invmap(n),candidates(n))
       allocate(t1(n),t2(n),t1tilde(n),t2tilde(n))
       call vf2pp_matching_order(a,b,order)
!
       map(:) = 0
       invmap(:) = 0
       t1(:) = .FALSE.
       t2(:) = .FALSE.
       t1tilde(:) = .TRUE.
       t2tilde(:) = .TRUE.
       found = .FALSE.
       call backtrack(1)
       iso = found
!
       deallocate(order,map,invmap,candidates,t1,t2,t1tilde,t2tilde)
!
       return
!
       contains
!
       recursive subroutine backtrack(level)
!
       integer,intent(in)  ::  level
       integer             ::  u,v,idx,ncand
!
       if ( found ) return
       if ( level .gt. n ) then
         found = .TRUE.
         return
       end if
!
       u = order(level)
       call vf2pp_find_candidates_undir(u,a,b,map,invmap,t2tilde,      &
                                        candidates,ncand)
       do idx = 1, ncand
         v = candidates(idx)
         if ( .not. vf2pp_feasible_pair_undir(u,v,a,b,map) ) cycle
         if ( .not. vf2pp_look_ahead_undir(u,v,a,b,t1,t2,t1tilde,     &
                                           t2tilde) ) cycle
         map(u) = v
         invmap(v) = u
         call vf2pp_update_tinout_undir(u,v,a,b,map,invmap,t1,t2,     &
                                        t1tilde,t2tilde)
         call backtrack(level+1)
         if ( found ) return
         map(u) = 0
         invmap(v) = 0
         call vf2pp_restore_tinout_undir(u,v,a,b,map,invmap,t1,t2,    &
                                         t1tilde,t2tilde)
       end do
!
       return
       end subroutine backtrack
!
       end function are_isomorphic_undir
!
!======================================================================!
!
       logical function are_isomorphic_dir(a,b) result(iso)
!
       type(iso_graph),intent(in)  ::  a,b
       integer                     ::  n
       integer,dimension(:),allocatable  ::  order,map,invmap
       integer,dimension(:),allocatable  ::  candidates
       logical,dimension(:),allocatable  ::  t1,t2,t1tilde,t2tilde
       logical,dimension(:),allocatable  ::  t1in,t2in
       logical                     ::  found
!
       n = a%n
       allocate(order(n),map(n),invmap(n),candidates(n))
       allocate(t1(n),t2(n),t1tilde(n),t2tilde(n),t1in(n),t2in(n))
       call vf2pp_matching_order(a,b,order)
!
       map(:) = 0
       invmap(:) = 0
       t1(:) = .FALSE.
       t2(:) = .FALSE.
       t1in(:) = .FALSE.
       t2in(:) = .FALSE.
       t1tilde(:) = .TRUE.
       t2tilde(:) = .TRUE.
       found = .FALSE.
       call backtrack(1)
       iso = found
!
       deallocate(order,map,invmap,candidates,t1,t2,t1tilde,t2tilde,  &
                  t1in,t2in)
!
       return
!
       contains
!
       recursive subroutine backtrack(level)
!
       integer,intent(in)  ::  level
       integer             ::  u,v,idx,ncand
!
       if ( found ) return
       if ( level .gt. n ) then
         found = .TRUE.
         return
       end if
!
       u = order(level)
       call vf2pp_find_candidates_dir(u,a,b,map,invmap,t2tilde,       &
                                      candidates,ncand)
       do idx = 1, ncand
         v = candidates(idx)
         if ( .not. vf2pp_feasible_pair_dir(u,v,a,b,map) ) cycle
         if ( .not. vf2pp_look_ahead_dir(u,v,a,b,t1,t2,t1in,t2in,     &
                                         t1tilde,t2tilde) ) cycle
         map(u) = v
         invmap(v) = u
         call vf2pp_update_tinout_dir(u,v,a,b,map,invmap,t1,t2,t1in,  &
                                      t2in,t1tilde,t2tilde)
         call backtrack(level+1)
         if ( found ) return
         map(u) = 0
         invmap(v) = 0
         call vf2pp_restore_tinout_dir(u,v,a,b,map,invmap,t1,t2,t1in, &
                                       t2in,t1tilde,t2tilde)
       end do
!
       return
       end subroutine backtrack
!
       end function are_isomorphic_dir
!
!======================================================================!
!
       subroutine vf2pp_matching_order(sg,fg,order)
!
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(out)   ::  order
       logical,dimension(:),allocatable   ::  unordered
       logical,dimension(:),allocatable   ::  inlayer
       logical,dimension(:),allocatable   ::  selected
       integer,dimension(:),allocatable   ::  used_degree
       integer,dimension(:),allocatable   ::  queue
       integer,dimension(:),allocatable   ::  level
       integer                            ::  n
       integer                            ::  norder
       integer                            ::  start
       integer                            ::  maxlevel
       integer                            ::  qhead,qtail
       integer                            ::  d,u,v,next
!
       n = sg%n
       allocate(unordered(n),inlayer(n),selected(n),used_degree(n),    &
                queue(n),level(n))
       unordered(:) = .TRUE.
       selected(:) = .FALSE.
       used_degree(:) = 0
       order(:) = 0
       norder = 0
!
       do while ( norder .lt. n )
         start = vf2pp_component_start(sg,fg,unordered,selected)
         inlayer(:) = .FALSE.
         level(:) = 0
         queue(:) = 0
         qhead = 1
         qtail = 1
         queue(qtail) = start
         level(start) = 1
         maxlevel = 1
!
         do while ( qhead .le. qtail )
           u = queue(qhead)
           qhead = qhead + 1
           do v = 1, n
             if ( .not. unordered(v) ) cycle
             if ( level(v) .ne. 0 ) cycle
             if ( vf2pp_connected(sg,u,v) ) then
               qtail = qtail + 1
               queue(qtail) = v
               level(v) = level(u) + 1
               if ( level(v) .gt. maxlevel ) maxlevel = level(v)
             end if
           end do
         end do
!
         do d = 1, maxlevel
           inlayer(:) = .FALSE.
           do u = 1, n
             if ( unordered(u) .and. level(u) .eq. d ) inlayer(u) =    &
                  .TRUE.
           end do
           do while ( any(inlayer) )
             next = vf2pp_best_layer_node(sg,fg,unordered,inlayer,     &
                                          selected,used_degree)
             if ( next .eq. 0 ) exit
             norder = norder + 1
             order(norder) = next
             unordered(next) = .FALSE.
             selected(next) = .TRUE.
             inlayer(next) = .FALSE.
             do v = 1, n
               if ( vf2pp_connected(sg,next,v) ) then
                 used_degree(v) = used_degree(v) + 1
               end if
             end do
           end do
         end do
       end do
!
       deallocate(unordered,inlayer,selected,used_degree,queue,level)
!
       return
       end subroutine vf2pp_matching_order
!
!======================================================================!
!
       integer function vf2pp_component_start(sg,fg,unordered,         &
                                              selected)                &
                         result(start)
!
       type(iso_graph),intent(in)       ::  sg,fg
       logical,dimension(:),intent(in)  ::  unordered
       logical,dimension(:),intent(in)  ::  selected
       integer                          ::  i
!
       start = 0
       do i = 1, sg%n
         if ( .not. unordered(i) ) cycle
         if ( start .eq. 0 ) then
           start = i
         else if ( vf2pp_better_start(sg,fg,i,start,selected) ) then
           start = i
         end if
       end do
!
       return
       end function vf2pp_component_start
!
!======================================================================!
!
       logical function vf2pp_better_start(sg,fg,u,v,selected)        &
                        result(better)
!
       type(iso_graph),intent(in)  ::  sg,fg
       integer,intent(in)          ::  u,v
       logical,dimension(:),intent(in) :: selected
       integer                     ::  ru,rv,du,dv
!
       ru = vf2pp_label_rarity(sg,fg,sg%labels(u),selected)
       rv = vf2pp_label_rarity(sg,fg,sg%labels(v),selected)
       du = vf2pp_order_degree(sg,u)
       dv = vf2pp_order_degree(sg,v)
!
       better = .FALSE.
       if ( ru .lt. rv ) then
         better = .TRUE.
       else if ( ru .eq. rv .and. du .gt. dv ) then
         better = .TRUE.
       end if
!
       return
       end function vf2pp_better_start
!
!======================================================================!
!
       integer function vf2pp_best_layer_node(sg,fg,unordered,inlayer, &
                                              selected,used_degree)    &
                         result(best)
!
       type(iso_graph),intent(in)       ::  sg,fg
       logical,dimension(:),intent(in)  ::  unordered
       logical,dimension(:),intent(in)  ::  inlayer
       logical,dimension(:),intent(in)  ::  selected
       integer,dimension(:),intent(in)  ::  used_degree
       integer                          ::  i
!
       best = 0
       do i = 1, sg%n
         if ( .not. unordered(i) .or. .not. inlayer(i) ) cycle
         if ( best .eq. 0 ) then
           best = i
         else if ( vf2pp_better_layer_choice(sg,fg,i,best,             &
                                             selected,used_degree) ) then
           best = i
         end if
       end do
!
       return
       end function vf2pp_best_layer_node
!
!======================================================================!
!
       logical function vf2pp_better_layer_choice(sg,fg,u,v,           &
                                                  selected,            &
                                                  used_degree)         &
                       result(better)
!
       type(iso_graph),intent(in)       ::  sg,fg
       integer,intent(in)               ::  u,v
       logical,dimension(:),intent(in)   ::  selected
       integer,dimension(:),intent(in)  ::  used_degree
       integer                          ::  ru,rv,du,dv
!
       ru = vf2pp_label_rarity(sg,fg,sg%labels(u),selected)
       rv = vf2pp_label_rarity(sg,fg,sg%labels(v),selected)
       du = vf2pp_order_degree(sg,u)
       dv = vf2pp_order_degree(sg,v)
!
       better = .FALSE.
       if ( used_degree(u) .gt. used_degree(v) ) then
         better = .TRUE.
       else if ( used_degree(u) .eq. used_degree(v) ) then
         if ( du .gt. dv ) then
           better = .TRUE.
         else if ( du .eq. dv .and. ru .lt. rv ) then
           better = .TRUE.
         end if
       end if
!
       return
       end function vf2pp_better_layer_choice
!
!======================================================================!
!
       integer function vf2pp_label_count(graph,label) result(nlabel)
!
       type(iso_graph),intent(in)  ::  graph
       character(len=*),intent(in) ::  label
       integer                     ::  i
!
       nlabel = 0
       do i = 1, graph%n
         if ( trim(graph%labels(i)) .eq. trim(label) ) nlabel =        &
              nlabel + 1
       end do
!
       return
       end function vf2pp_label_count
!
!======================================================================!
!
       integer function vf2pp_label_rarity(sg,fg,label,selected)       &
                         result(rarity)
!
       type(iso_graph),intent(in)       ::  sg,fg
       character(len=*),intent(in)      ::  label
       logical,dimension(:),intent(in)  ::  selected
       integer                          ::  i
!
       rarity = vf2pp_label_count(fg,label)
       do i = 1, sg%n
         if ( selected(i) ) then
           if ( trim(sg%labels(i)) .eq. trim(label) ) rarity =         &
                rarity - 1
         end if
       end do
!
       return
       end function vf2pp_label_rarity
!
!======================================================================!
!
       logical function vf2pp_connected(graph,i,j) result(connected)
!
       type(iso_graph),intent(in)  ::  graph
       integer,intent(in)          ::  i,j
!
       connected = graph%adj(i,j)
       if ( graph%directed ) connected = connected .or. graph%adj(j,i)
!
       return
       end function vf2pp_connected
!
!======================================================================!
!
       integer function vf2pp_order_degree(graph,node) result(odeg)
!
       type(iso_graph),intent(in)  ::  graph
       integer,intent(in)          ::  node
       integer                     ::  i
!
       if ( .not. graph%directed ) then
         odeg = graph%degree(node)
         return
       end if
!
       odeg = 0
       do i = 1, graph%n
         if ( i .eq. node ) cycle
         if ( graph%adj(node,i) .or. graph%adj(i,node) ) odeg =       &
              odeg + 1
       end do
!
       return
       end function vf2pp_order_degree
!
!======================================================================!
!
       subroutine vf2pp_find_candidates_undir(u,sg,fg,map,invmap,      &
                                              t2tilde,candidates,      &
                                              ncand)
!
       integer,intent(in)                 ::  u
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(in)    ::  t2tilde
       integer,dimension(:),intent(out)   ::  candidates
       integer,intent(out)                ::  ncand
       integer                            ::  v,w
       logical,dimension(:),allocatable   ::  candmask
       logical                            ::  has_covered
!
       ncand = 0
       has_covered = .FALSE.
       allocate(candmask(fg%n))
       candmask(:) = .FALSE.
       do w = 1, sg%n
         if ( map(w) .eq. 0 ) cycle
         if ( .not. sg%adj(u,w) ) cycle
         if ( .not. has_covered ) then
           candmask(:) = fg%adj(map(w),:)
           has_covered = .TRUE.
         else
           candmask(:) = candmask(:) .and. fg%adj(map(w),:)
         end if
       end do
       if ( .not. has_covered ) candmask(:) = t2tilde(:)
!
       do v = 1, fg%n
         if ( .not. candmask(v) ) cycle
         if ( invmap(v) .ne. 0 ) cycle
         if ( .not. same_labeled_vertex_invariant_undir(sg,u,fg,v) )   &
              cycle
         ncand = ncand + 1
         candidates(ncand) = v
       end do
       deallocate(candmask)
!
       return
       end subroutine vf2pp_find_candidates_undir
!
!======================================================================!
!
       subroutine vf2pp_find_candidates_dir(u,sg,fg,map,invmap,        &
                                            t2tilde,candidates,ncand)
!
       integer,intent(in)                 ::  u
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(in)    ::  t2tilde
       integer,dimension(:),intent(out)   ::  candidates
       integer,intent(out)                ::  ncand
       integer                            ::  v,w
       logical,dimension(:),allocatable   ::  candmask
       logical,dimension(:),allocatable   ::  tmpmask
       logical                            ::  has_covered
!
       ncand = 0
       has_covered = .FALSE.
       allocate(candmask(fg%n),tmpmask(fg%n))
       candmask(:) = .FALSE.
       do w = 1, sg%n
         if ( map(w) .eq. 0 ) cycle
         if ( sg%adj(u,w) ) then
           tmpmask(:) = fg%adj(:,map(w))
           if ( .not. has_covered ) then
             candmask(:) = tmpmask(:)
             has_covered = .TRUE.
           else
             candmask(:) = candmask(:) .and. tmpmask(:)
           end if
         end if
         if ( sg%adj(w,u) ) then
           tmpmask(:) = fg%adj(map(w),:)
           if ( .not. has_covered ) then
             candmask(:) = tmpmask(:)
             has_covered = .TRUE.
           else
             candmask(:) = candmask(:) .and. tmpmask(:)
           end if
         end if
       end do
       if ( .not. has_covered ) candmask(:) = t2tilde(:)
!
       do v = 1, fg%n
         if ( .not. candmask(v) ) cycle
         if ( invmap(v) .ne. 0 ) cycle
         if ( .not. same_labeled_vertex_invariant_dir(sg,u,fg,v) )     &
              cycle
         ncand = ncand + 1
         candidates(ncand) = v
       end do
       deallocate(candmask,tmpmask)
!
       return
       end subroutine vf2pp_find_candidates_dir
!
!======================================================================!
!
       logical function vf2pp_feasible_pair_undir(u,v,sg,fg,map)       &
                        result(feasible)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map
       integer                            ::  w
!
       feasible = .TRUE.
       do w = 1, sg%n
         if ( map(w) .eq. 0 ) cycle
         if ( sg%adj(u,w) .neqv. fg%adj(v,map(w)) ) then
           feasible = .FALSE.
           return
         end if
       end do
!
       return
       end function vf2pp_feasible_pair_undir
!
!======================================================================!
!
       logical function vf2pp_feasible_pair_dir(u,v,sg,fg,map)         &
                        result(feasible)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map
       integer                            ::  w
!
       feasible = .TRUE.
       do w = 1, sg%n
         if ( map(w) .eq. 0 ) cycle
         if ( sg%adj(u,w) .neqv. fg%adj(v,map(w)) ) then
           feasible = .FALSE.
           return
         end if
         if ( sg%adj(w,u) .neqv. fg%adj(map(w),v) ) then
           feasible = .FALSE.
           return
         end if
       end do
!
       return
       end function vf2pp_feasible_pair_dir
!
!======================================================================!
!
       logical function vf2pp_look_ahead_undir(u,v,sg,fg,t1,t2,        &
                                               t1tilde,t2tilde)        &
                        result(feasible)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       logical,dimension(:),intent(in)    ::  t1,t2,t1tilde,t2tilde
!
       feasible = vf2pp_labcut_undir_part(u,v,sg,fg,t1,t2)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_undir_part(u,v,sg,fg,t1tilde,t2tilde)
!
       return
       end function vf2pp_look_ahead_undir
!
!======================================================================!
!
       logical function vf2pp_look_ahead_dir(u,v,sg,fg,t1,t2,t1in,     &
                                             t2in,t1tilde,t2tilde)     &
                        result(feasible)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       logical,dimension(:),intent(in)    ::  t1,t2,t1in,t2in
       logical,dimension(:),intent(in)    ::  t1tilde,t2tilde
!
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1,t2,.TRUE.)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1in,t2in,.TRUE.)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1tilde,t2tilde,    &
                                        .TRUE.)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1,t2,.FALSE.)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1in,t2in,.FALSE.)
       if ( .not. feasible ) return
       feasible = vf2pp_labcut_dir_part(u,v,sg,fg,t1tilde,t2tilde,    &
                                        .FALSE.)
!
       return
       end function vf2pp_look_ahead_dir
!
!======================================================================!
!
       logical function vf2pp_labcut_undir_part(u,v,sg,fg,front1,     &
                                                front2) result(ok)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       logical,dimension(:),intent(in)    ::  front1,front2
       character(len=leniso),dimension(:),allocatable ::  labels
       integer,dimension(:),allocatable   ::  counts
       integer                            ::  nlabels
       integer                            ::  i
!
       allocate(labels(sg%n+fg%n),counts(sg%n+fg%n))
       labels(:) = ''
       counts(:) = 0
       nlabels = 0
!
       do i = 1, sg%n
         if ( front1(i) .and. sg%adj(u,i) ) call vf2pp_add_label_delta(&
              labels,counts,nlabels,sg%labels(i),1)
       end do
       do i = 1, fg%n
         if ( front2(i) .and. fg%adj(v,i) ) call vf2pp_add_label_delta(&
              labels,counts,nlabels,fg%labels(i),-1)
       end do
!
       ok = vf2pp_zero_label_deltas(counts,nlabels)
       deallocate(labels,counts)
!
       return
       end function vf2pp_labcut_undir_part
!
!======================================================================!
!
       logical function vf2pp_labcut_dir_part(u,v,sg,fg,front1,front2,&
                                              forward) result(ok)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       logical,dimension(:),intent(in)    ::  front1,front2
       logical,intent(in)                 ::  forward
       character(len=leniso),dimension(:),allocatable ::  labels
       integer,dimension(:),allocatable   ::  counts
       integer                            ::  nlabels
       integer                            ::  i
!
       allocate(labels(sg%n+fg%n),counts(sg%n+fg%n))
       labels(:) = ''
       counts(:) = 0
       nlabels = 0
!
       do i = 1, sg%n
         if ( .not. front1(i) ) cycle
         if ( forward ) then
           if ( sg%adj(u,i) ) call vf2pp_add_label_delta(labels,       &
                counts,nlabels,sg%labels(i),1)
         else
           if ( sg%adj(i,u) ) call vf2pp_add_label_delta(labels,       &
                counts,nlabels,sg%labels(i),1)
         end if
       end do
       do i = 1, fg%n
         if ( .not. front2(i) ) cycle
         if ( forward ) then
           if ( fg%adj(v,i) ) call vf2pp_add_label_delta(labels,       &
                counts,nlabels,fg%labels(i),-1)
         else
           if ( fg%adj(i,v) ) call vf2pp_add_label_delta(labels,       &
                counts,nlabels,fg%labels(i),-1)
         end if
       end do
!
       ok = vf2pp_zero_label_deltas(counts,nlabels)
       deallocate(labels,counts)
!
       return
       end function vf2pp_labcut_dir_part
!
!======================================================================!
!
       subroutine vf2pp_add_label_delta(labels,counts,nlabels,label,   &
                                        delta)
!
       character(len=leniso),dimension(:),intent(inout) ::  labels
       integer,dimension(:),intent(inout)               ::  counts
       integer,intent(inout)                            ::  nlabels
       character(len=*),intent(in)                      ::  label
       integer,intent(in)                               ::  delta
       integer                                          ::  i
!
       do i = 1, nlabels
         if ( trim(labels(i)) .eq. trim(label) ) then
           counts(i) = counts(i) + delta
           return
         end if
       end do
!
       nlabels = nlabels + 1
       labels(nlabels) = label
       counts(nlabels) = delta
!
       return
       end subroutine vf2pp_add_label_delta
!
!======================================================================!
!
       logical function vf2pp_zero_label_deltas(counts,nlabels)        &
                        result(ok)
!
       integer,dimension(:),intent(in)  ::  counts
       integer,intent(in)               ::  nlabels
       integer                          ::  i
!
       ok = .TRUE.
       do i = 1, nlabels
         if ( counts(i) .ne. 0 ) then
           ok = .FALSE.
           return
         end if
       end do
!
       return
       end function vf2pp_zero_label_deltas
!
!======================================================================!
!
       subroutine vf2pp_update_tinout_undir(u,v,sg,fg,map,invmap,t1,  &
                                            t2,t1tilde,t2tilde)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(inout) ::  t1,t2,t1tilde,t2tilde
       integer                            ::  i
!
       do i = 1, sg%n
         if ( .not. sg%adj(u,i) ) cycle
         if ( map(i) .ne. 0 ) cycle
         t1(i) = .TRUE.
         t1tilde(i) = .FALSE.
       end do
       do i = 1, fg%n
         if ( .not. fg%adj(v,i) ) cycle
         if ( invmap(i) .ne. 0 ) cycle
         t2(i) = .TRUE.
         t2tilde(i) = .FALSE.
       end do
!
       t1(u) = .FALSE.
       t2(v) = .FALSE.
       t1tilde(u) = .FALSE.
       t2tilde(v) = .FALSE.
!
       return
       end subroutine vf2pp_update_tinout_undir
!
!======================================================================!
!
       subroutine vf2pp_restore_tinout_undir(u,v,sg,fg,map,invmap,t1, &
                                             t2,t1tilde,t2tilde)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(inout) ::  t1,t2,t1tilde,t2tilde
       integer                            ::  i
       logical                            ::  added
!
       added = .FALSE.
       do i = 1, sg%n
         if ( .not. sg%adj(u,i) ) cycle
         if ( map(i) .ne. 0 ) then
           added = .TRUE.
           t1(u) = .TRUE.
           t1tilde(u) = .FALSE.
         else if ( .not. vf2pp_has_mapped_neighbor_undir(sg,i,map) ) then
           t1(i) = .FALSE.
           t1tilde(i) = .TRUE.
         end if
       end do
       if ( .not. added ) then
         t1(u) = .FALSE.
         t1tilde(u) = .TRUE.
       end if
!
       added = .FALSE.
       do i = 1, fg%n
         if ( .not. fg%adj(v,i) ) cycle
         if ( invmap(i) .ne. 0 ) then
           added = .TRUE.
           t2(v) = .TRUE.
           t2tilde(v) = .FALSE.
         else if ( .not. vf2pp_has_mapped_neighbor_undir(fg,i,        &
                    invmap) ) then
           t2(i) = .FALSE.
           t2tilde(i) = .TRUE.
         end if
       end do
       if ( .not. added ) then
         t2(v) = .FALSE.
         t2tilde(v) = .TRUE.
       end if
!
       return
       end subroutine vf2pp_restore_tinout_undir
!
!======================================================================!
!
       subroutine vf2pp_update_tinout_dir(u,v,sg,fg,map,invmap,t1,t2, &
                                          t1in,t2in,t1tilde,t2tilde)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(inout) ::  t1,t2,t1in,t2in
       logical,dimension(:),intent(inout) ::  t1tilde,t2tilde
       integer                            ::  i
!
       do i = 1, sg%n
         if ( sg%adj(u,i) .and. map(i) .eq. 0 ) then
           t1(i) = .TRUE.
           t1tilde(i) = .FALSE.
         end if
         if ( sg%adj(i,u) .and. map(i) .eq. 0 ) then
           t1in(i) = .TRUE.
           t1tilde(i) = .FALSE.
         end if
       end do
       do i = 1, fg%n
         if ( fg%adj(v,i) .and. invmap(i) .eq. 0 ) then
           t2(i) = .TRUE.
           t2tilde(i) = .FALSE.
         end if
         if ( fg%adj(i,v) .and. invmap(i) .eq. 0 ) then
           t2in(i) = .TRUE.
           t2tilde(i) = .FALSE.
         end if
       end do
!
       t1(u) = .FALSE.
       t2(v) = .FALSE.
       t1in(u) = .FALSE.
       t2in(v) = .FALSE.
       t1tilde(u) = .FALSE.
       t2tilde(v) = .FALSE.
!
       return
       end subroutine vf2pp_update_tinout_dir
!
!======================================================================!
!
       subroutine vf2pp_restore_tinout_dir(u,v,sg,fg,map,invmap,t1,   &
                                           t2,t1in,t2in,t1tilde,      &
                                           t2tilde)
!
       integer,intent(in)                 ::  u,v
       type(iso_graph),intent(in)         ::  sg,fg
       integer,dimension(:),intent(in)    ::  map,invmap
       logical,dimension(:),intent(inout) ::  t1,t2,t1in,t2in
       logical,dimension(:),intent(inout) ::  t1tilde,t2tilde
       integer                            ::  i
       logical                            ::  added
!
       added = .FALSE.
       do i = 1, sg%n
         if ( sg%adj(u,i) ) then
           if ( map(i) .ne. 0 ) then
             added = .TRUE.
             t1in(u) = .TRUE.
             t1tilde(u) = .FALSE.
           else
             call vf2pp_restore_dir_node(sg,i,map,t1,t1in,t1tilde)
           end if
         end if
         if ( sg%adj(i,u) ) then
           if ( map(i) .ne. 0 ) then
             added = .TRUE.
             t1(u) = .TRUE.
             t1tilde(u) = .FALSE.
           else
             call vf2pp_restore_dir_node(sg,i,map,t1,t1in,t1tilde)
           end if
         end if
       end do
       if ( .not. added ) then
         t1(u) = .FALSE.
         t1in(u) = .FALSE.
         t1tilde(u) = .TRUE.
       end if
!
       added = .FALSE.
       do i = 1, fg%n
         if ( fg%adj(v,i) ) then
           if ( invmap(i) .ne. 0 ) then
             added = .TRUE.
             t2in(v) = .TRUE.
             t2tilde(v) = .FALSE.
           else
             call vf2pp_restore_dir_node(fg,i,invmap,t2,t2in,         &
                                         t2tilde)
           end if
         end if
         if ( fg%adj(i,v) ) then
           if ( invmap(i) .ne. 0 ) then
             added = .TRUE.
             t2(v) = .TRUE.
             t2tilde(v) = .FALSE.
           else
             call vf2pp_restore_dir_node(fg,i,invmap,t2,t2in,         &
                                         t2tilde)
           end if
         end if
       end do
       if ( .not. added ) then
         t2(v) = .FALSE.
         t2in(v) = .FALSE.
         t2tilde(v) = .TRUE.
       end if
!
       return
       end subroutine vf2pp_restore_tinout_dir
!
!======================================================================!
!
       subroutine vf2pp_restore_dir_node(graph,node,mapped,tout,tin,   &
                                         ttilde)
!
       type(iso_graph),intent(in)         ::  graph
       integer,intent(in)                 ::  node
       integer,dimension(:),intent(in)    ::  mapped
       logical,dimension(:),intent(inout) ::  tout,tin,ttilde
!
       if ( .not. vf2pp_has_mapped_predecessor_dir(graph,node,        &
            mapped) ) tout(node) = .FALSE.
       if ( .not. vf2pp_has_mapped_successor_dir(graph,node,mapped) ) &
            tin(node) = .FALSE.
       if ( .not. tout(node) .and. .not. tin(node) ) ttilde(node) =   &
            .TRUE.
!
       return
       end subroutine vf2pp_restore_dir_node
!
!======================================================================!
!
       logical function vf2pp_has_mapped_neighbor_undir(graph,node,    &
                                                        mapped)        &
                        result(has_mapped)
!
       type(iso_graph),intent(in)       ::  graph
       integer,intent(in)               ::  node
       integer,dimension(:),intent(in)  ::  mapped
       integer                          ::  i
!
       has_mapped = .FALSE.
       do i = 1, graph%n
         if ( graph%adj(node,i) .and. mapped(i) .ne. 0 ) then
           has_mapped = .TRUE.
           return
         end if
       end do
!
       return
       end function vf2pp_has_mapped_neighbor_undir
!
!======================================================================!
!
       logical function vf2pp_has_mapped_successor_dir(graph,node,     &
                                                       mapped)         &
                        result(has_mapped)
!
       type(iso_graph),intent(in)       ::  graph
       integer,intent(in)               ::  node
       integer,dimension(:),intent(in)  ::  mapped
       integer                          ::  i
!
       has_mapped = .FALSE.
       do i = 1, graph%n
         if ( graph%adj(node,i) .and. mapped(i) .ne. 0 ) then
           has_mapped = .TRUE.
           return
         end if
       end do
!
       return
       end function vf2pp_has_mapped_successor_dir
!
!======================================================================!
!
       logical function vf2pp_has_mapped_predecessor_dir(graph,node,   &
                                                         mapped)       &
                        result(has_mapped)
!
       type(iso_graph),intent(in)       ::  graph
       integer,intent(in)               ::  node
       integer,dimension(:),intent(in)  ::  mapped
       integer                          ::  i
!
       has_mapped = .FALSE.
       do i = 1, graph%n
         if ( graph%adj(i,node) .and. mapped(i) .ne. 0 ) then
           has_mapped = .TRUE.
           return
         end if
       end do
!
       return
       end function vf2pp_has_mapped_predecessor_dir
!
!======================================================================!
!
       logical function same_tag_multiset(a,b) result(same)
!
       type(iso_graph),intent(in)  ::  a,b
       logical,dimension(:),allocatable ::  used
       integer                     ::  i,j,n
!
       n = a%n
       allocate(used(n))
       used(:) = .FALSE.
       same = .TRUE.
       do i = 1, n
         do j = 1, n
           if ( .not. used(j) ) then
             if ( same_labeled_vertex_invariant(a,i,b,j) ) then
               used(j) = .TRUE.
               exit
             end if
           end if
         end do
         if ( j .gt. n ) then
           same = .FALSE.
           exit
         end if
       end do
       deallocate(used)
!
       return
       end function same_tag_multiset
!
!======================================================================!
!
       logical function same_labeled_vertex_invariant(a,ia,b,ib)       &
                        result(same)
!
       type(iso_graph),intent(in)  ::  a,b
       integer,intent(in)          ::  ia,ib
!
       same = .FALSE.
!      Label preservation is the defining constraint of the aggregate
!      labelled-graph isomorphism: a body/group can only be mapped onto
!      another body/group with the same name.
       if ( trim(a%labels(ia)) .ne. trim(b%labels(ib)) ) return
!      Degree checks are graph invariants used to prune impossible labelled
!      mappings before entering the VF2++ search tree.
       if ( a%directed ) then
         if ( a%indegree(ia) .ne. b%indegree(ib) ) return
         if ( a%outdegree(ia) .ne. b%outdegree(ib) ) return
       else
         if ( a%degree(ia) .ne. b%degree(ib) ) return
       end if
       same = .TRUE.
!
       return
       end function same_labeled_vertex_invariant
!
!======================================================================!
!
       logical function same_labeled_vertex_invariant_undir(a,ia,b,ib) &
                        result(same)
!
       type(iso_graph),intent(in)  ::  a,b
       integer,intent(in)          ::  ia,ib
!
       same = .FALSE.
       if ( trim(a%labels(ia)) .ne. trim(b%labels(ib)) ) return
       if ( a%degree(ia) .ne. b%degree(ib) ) return
       same = .TRUE.
!
       return
       end function same_labeled_vertex_invariant_undir
!
!======================================================================!
!
       logical function same_labeled_vertex_invariant_dir(a,ia,b,ib)   &
                        result(same)
!
       type(iso_graph),intent(in)  ::  a,b
       integer,intent(in)          ::  ia,ib
!
       same = .FALSE.
       if ( trim(a%labels(ia)) .ne. trim(b%labels(ib)) ) return
       if ( a%indegree(ia) .ne. b%indegree(ib) ) return
       if ( a%outdegree(ia) .ne. b%outdegree(ib) ) return
       same = .TRUE.
!
       return
       end function same_labeled_vertex_invariant_dir
!
!======================================================================!
!
       subroutine add_iso_class(bucket,graph,class_id)
!
       type(iso_bucket),intent(inout)     ::  bucket
       type(iso_graph),intent(in)         ::  graph
       integer,intent(out)                ::  class_id
       integer                            ::  n
!
       n = bucket%nclasses
       if ( n .eq. bucket%class_capacity ) call grow_iso_classes(bucket)
       bucket%classes(n+1)%representative = graph
       bucket%classes(n+1)%count = 1
       bucket%classes(n+1)%old_id = n + 1
       bucket%nclasses = n + 1
       class_id = n + 1
!
       return
       end subroutine add_iso_class
!
!======================================================================!
!
       subroutine grow_iso_classes(bucket)
!
       type(iso_bucket),intent(inout)              ::  bucket
       type(iso_class),dimension(:),allocatable    ::  tmp
       integer                                     ::  i
       integer                                     ::  new_capacity
!
       if ( bucket%class_capacity .ge. max_iso_class_capacity ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Maximum number of isomorphism '// &
                           'classes exceeded'
         write(*,*)
         write(*,'(3X,A,I0)') 'Aggregate identifier : ',bucket%iagg
         write(*,'(3X,A,I0)') 'Maximum classes      : ',               &
                              max_iso_class_capacity
         write(*,*)
         write(*,'(3X,A)') 'Increase --max-iso-classes or revise '//  &
                           'the connectivity criteria.'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       new_capacity = min(2*bucket%class_capacity,                     &
                          max_iso_class_capacity)
       allocate(tmp(new_capacity))
       do i = 1, bucket%nclasses
         tmp(i) = bucket%classes(i)
       end do
       call move_alloc(tmp,bucket%classes)
       bucket%class_capacity = new_capacity
!
       return
       end subroutine grow_iso_classes
!
!======================================================================!
!
       subroutine append_assignment(bucket,class_id)
!
       type(iso_bucket),intent(inout)     ::  bucket
       integer,intent(in)                 ::  class_id
!
       if ( bucket%nobs .gt. bucket%assignment_capacity )              &
         call grow_iso_assignments(bucket)
       bucket%assignments(bucket%nobs) = class_id
!
       return
       end subroutine append_assignment
!
!======================================================================!
!
       subroutine grow_iso_assignments(bucket)
!
       type(iso_bucket),intent(inout)             ::  bucket
       integer,dimension(:),allocatable           ::  tmp
       integer                                    ::  i
       integer                                    ::  new_capacity
!
       new_capacity = max(2*bucket%assignment_capacity,bucket%nobs)
       allocate(tmp(new_capacity))
       do i = 1, bucket%nobs-1
         tmp(i) = bucket%assignments(i)
       end do
       call move_alloc(tmp,bucket%assignments)
       bucket%assignment_capacity = new_capacity
!
       return
       end subroutine grow_iso_assignments
!
!======================================================================!
!
       subroutine write_iso_bucket(bucket)
!
       type(iso_bucket),intent(inout)  ::  bucket
       integer,dimension(:),allocatable ::  order
       integer,dimension(:),allocatable ::  old2new
       integer                         ::  i,j
       real(kind=8)                    ::  pop
       character(len=lenout)           ::  fname
!
       if ( bucket%nobs .le. 0 ) return
!
       allocate(order(bucket%nclasses),old2new(bucket%nclasses))
       do i = 1, bucket%nclasses
         order(i) = i
       end do
       call sort_class_order(bucket,order)
       do i = 1, bucket%nclasses
         old2new(order(i)) = i
       end do
!
       fname = trim(bucket%filebase)//'_iso_classes.dat'
       open(unit=91,file=trim(fname),action='write')
       write(91,'(A)') '# class_id count population_percent'
       do i = 1, bucket%nclasses
         j = order(i)
         pop = 100.0d0*dble(bucket%classes(j)%count)/dble(bucket%nobs)
         write(91,'(I8,1X,I12,1X,F16.8)') i,bucket%classes(j)%count,pop
       end do
       close(91)
!
       fname = trim(bucket%filebase)//'_iso_assignments.dat'
       open(unit=91,file=trim(fname),action='write')
       write(91,'(A)') '# observation_id class_id'
       do i = 1, bucket%nobs
         write(91,'(I12,1X,I8)') i,old2new(bucket%assignments(i))
       end do
       close(91)
!
       fname = trim(bucket%filebase)//'_iso_representatives.txt'
       open(unit=91,file=trim(fname),action='write')
       write(91,'(A)') trim(bucket%labels)
       do i = 1, bucket%nclasses
         j = order(i)
         pop = 100.0d0*dble(bucket%classes(j)%count)/dble(bucket%nobs)
         write(91,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,F16.8)')               &
              '# class_id',i,'count',bucket%classes(j)%count,          &
              'population_percent',pop
         call write_graph_matrix(91,bucket%classes(j)%representative)
       end do
       close(91)
!
       deallocate(order,old2new)
!
       return
      end subroutine write_iso_bucket
!
!======================================================================!
!
       subroutine iso_filebase(filebase,base)
!
       character(len=*),intent(in)      ::  filebase
       character(len=lenout),intent(out) ::  base
       integer                          ::  n
!
       base = filebase
       n = len_trim(base)
       if ( n .ge. 4 ) then
         if ( base(n-3:n) .eq. '.txt' ) base = base(:n-4)
       end if
!
       return
       end subroutine iso_filebase
!
!======================================================================!
!
       subroutine sort_class_order(bucket,order)
!
       type(iso_bucket),intent(in)       ::  bucket
       integer,dimension(:),intent(inout) ::  order
       integer                           ::  i,j,best,tmp
!
       do i = 1, size(order)-1
         best = i
         do j = i+1, size(order)
           if ( bucket%classes(order(j))%count .gt.                    &
                bucket%classes(order(best))%count ) best = j
         end do
         if ( best .ne. i ) then
           tmp = order(i)
           order(i) = order(best)
           order(best) = tmp
         end if
       end do
!
       return
       end subroutine sort_class_order
!
!======================================================================!
!
       subroutine write_graph_matrix(unit,graph)
!
       integer,intent(in)          ::  unit
       type(iso_graph),intent(in)  ::  graph
       integer                     ::  i,j
!
       do i = 1, graph%n
         write(unit,*) (graph%adj(i,j),j=1,graph%n)
       end do
!
       return
       end subroutine write_graph_matrix
!
       end module isotools
!
!======================================================================!
