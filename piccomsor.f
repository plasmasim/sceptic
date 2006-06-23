c Here are stored the common declarations for the parallel bloc solver



c     Iteration counter and maximum nb of iterations
      integer isor_k,isor_mi

c     Iteration parameters
      real sor_jac,sor_eps,sor_del

      common /sor_decl/ isor_k,isor_mi,sor_jac,sor_eps,sor_del


cccc    Parallel data    ccccc
c ---------------------------------------------

c     Number of participating processors
      integer nprocs

c     Process ID in cartesian topology
      integer mycartid

c     Cartesian communicator
      integer icommcart

c     dimension of each topology (nb of blocks)
      integer idims(2)

c     cartesian topology coords of this process (block)
      integer icoords(2)

c     origin of the blocks (in the total i,j frame), in a linear referencing
      integer iorig(5,5)

c     side length of each block (in the total i,j frame) for each dimension
c     iside(*,1) is the general value, iside(*,2) the uppermost value
      integer iside(2,2)

c     dimension of my block (in the total i,j frame)
      integer myside(2)

c     Whether the topology dimension is periodic or not
      logical lperiod(2)

c     iface is a new datatype containing the indexes of the points
c     to transmit to the neighbors on the linear array u
c     Arguments : ndims, parity, bulk or top
      integer iface(2,2,2)

      integer ktype(20)

      integer ibt(2)

      integer isdl(2),iddl(2),isdr(2),iddr(2)

      common /sor_parallel/ nprocs,icommcart,idims,icoords,iorig,
     $     iside,myside,lperiod,mycartid,iface,ktype,
     $     isdl,iddl,isdr,iddr,ibt

c ---------------------------------------------

c     Arrays for constructing Alltoall calls
      parameter (maxprocs=100)
      integer isdispls(0:maxprocs),irdispls(0:maxprocs)
      integer istypes(0:maxprocs),irtypes(0:maxprocs)
      integer iscounts(0:maxprocs),ircounts(0:maxprocs)
      integer iconp(10)
      integer isizeofreal,ith0

      common /alltoall/isdispls,irdispls,istypes,irtypes,
     $     iscounts,ircounts,iconp,isizeofreal,ith0
