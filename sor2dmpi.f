
c***************************************************
      subroutine sorparinit(myid2,sor_comm)

c     Creates a new communicator, sor_comm, who contains a subset of
c     MPI_COMM_WORLD, because we don't need all the processes for the
c     bloc SOR

c     Not to be confused with mpicommcart, which is a reorganisation of
c     sor_comm

c     myid2 is the process id in the new communicator. If myid2 is negative,
c     the process does not belong to the communicator

      integer myid2,nprocsor
      integer members(0:100)
      integer group_world,sor_group,sor_comm

      include 'piccom.f'
      include 'mpif.h'

      idim1=nint(sqrt(numprocs*NRUSED/NTHUSED*1.))
      idim2=nint(numprocs/idim1*1.)

      if((idim1.ge.numprocs)) then
         idim1=numprocs
         idim2=1
      elseif((idim2.ge.numprocs)) then
         idim2=numprocs
         idim1=1
      endif



      nprocsor=idim1*idim2

      do i=0,nprocsor-1
         members(i)=i
      enddo

      call MPI_COMM_GROUP(MPI_COMM_WORLD,group_world,ierr)
      call MPI_GROUP_INCL(group_world,nprocsor,members,sor_group,ierr)
      call MPI_GROUP_RANK(sor_group,myid2,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,sor_group,sor_comm,ierr)
      if (myid2.ge.0) then
         call MPI_COMM_RANK(sor_comm,myid2,ierr)
      endif

      if (myid2.eq.0) then
c         write(*,*) "sorp : ",idim1,idim2,nprocsor,numprocs
         write(*,510) idim1,idim2,nprocsor
 510     format('BlocX=',i2,'  BlocY=',i2,'  SorProc=',i2)
      endif

      end

c ------------------------------------------------------------------

c     Solve an elliptical problem in 2d by simultaneous overrelaxation
c     L(u) + f(u) = q(x,y), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c Ian Hutchinson, January 2006.
c
      subroutine sor2dmpi(sor_comm,Li,Lj,ni,nj,debyelen,bcphi,u,q,ictl,
     $     ierr,mpiid,idim1,idim2,apc,bpc,cpc,dpc,fpc,gpc)

      integer sor_comm,mpiid
c     The number of dimensions, 2 here.
      integer nd
      parameter (nd=2,nd2=nd*2)
c     Li leading dimension of i 
c     ni active dimension of i 
c     nj active dimension of j 
      integer Li,Lj,ni,nj
c     cij coefficients of the finite difference scheme, in the order
c         east,west,north,south [regarding i,j, as x,y].
c      real cij(nd2,Li,nj)
c     In this mpi routine we must use linear addressing for cij,u,q
c     so that the pointers can be used for each block.
      real apc(*),bpc(*),cpc(*),dpc(*),fpc(*)
      real gpc(0:Lj-1,1:5)
c     u  potential to be solved for (initialized on entry).
c        its boundaries are at 1,ni; 1,nj; 

c      real u(Li,nj)
      real u(*)
c     q  "charge density" input
c      real q(Li,nj)
      real q(*)
c User supplied functions which should be declared external in the 
c calling routine.

c     mpiid Returns my MPI process number for this mpi version.
      integer ictl,ierr
c     idim1,2   The number of blocks in dimensions 1 and 2.
      integer idim1,idim2
c Other things that we might want control over include the maximum
c number of iterations, the Jacobi radius, and the convergence size.


      common /sor2dctl/isor_mi,sor_jac,sor_eps,sor_del,isor_k

      real delta
      logical lconverged

c     Debyelength
      real debyelen

c     BC
      integer bcphi
c     lflag, decide if we call bbdy for the first time or not
      logical lflag
      
      include 'piccomsor.f'

c     Whether the topology dimension is periodic or not
      logical lperiod(ndims)
      data lperiod/ndims*.false./

c The origin of blocks structure may be considered
c      integer iorig(idim1+1,idim2+1)
c we declare it as a 1-d array. This is the first time here:
c It must be of dimension greater than the number of processes (blocks)
c      parameter (norigmax=1000)
c      integer iorig(norigmax)
c      parameter (ndims=2)
c bbdydecl declares most things for bbdy, using parameter ndims.
      
c-------------------------------------------------------------------
      
      save lflag

      idims(1)=idim1
      idims(2)=idim2
      iuds(1)=ni
      iuds(2)=nj
      ifull(1)=Li
      ifull(2)=nj
      if((idims(1)+1)*(idims(2)+1).gt.norigmax)then
         write(*,*)'Too many processes',idims(1),'x',idims(2),
     $        ' for norigmax=',norigmax
         stop
      endif
c Define mpi block structure.
      if(.not.lflag)  
     $     call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)

c End of control functions.
c-------------------------------------------------------------------
      ierr=0
      omega=1.
      lconverged=.false.
c     Main iteration

      do isor_k=1,isor_mi
         delta=0.
         umin=1.e20
         umax=0.
         relax=(omega*debyelen**2+1)/(debyelen**2+1)

c Do block boundary communications, returns block info icoords...myid.
c         if (mod(isor_k+2,3).eq.0) then
         call bbdy(sor_comm,iLs,iuds,u,isor_k,iorig,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,myorig1,myorig2,icommcart
     $        ,mpiid,lflag,out)


c         endif
c Do a relaxation.
         call sorrelax(isor_k,Li,Lj,myside(1),myside(2),
     $        u(myorig),q(myorig),apc(myorig1),bpc(myorig1),
     $        cpc(myorig),dpc(myorig),fpc(myorig),gpc(myorig2,1),out,
     $        relax,delta,bcphi)
         

c     Test convergence
         if (mod(isor_k+2,30).eq.0) then
            call testifconverged(sor_eps,delta, lconverged,sor_comm)
         endif
    
         if(lconverged.and.isor_k.ge.2)goto 11
c Chebychev acceleration:
         if(isor_k.eq.1)then
            omega=1./(1.-0.5*sor_jac**2)
         else
            omega=1./(1.-0.25*sor_jac**2*omega)
         endif
      enddo
      
c We finished the loop, implies we did not converge.
c      isor_k=-isor_mi

c-------------------------------------------------------------------
 11   continue

c Do the final mpi_gather [or allgather if all processes need
c the result].
      nk=-1
      

      call bbdy(sor_comm,iLs,iuds,u,nk,iorig,ndims,idims,lperiod,
     $     icoords,iLcoords,myside,myorig,myorig1,myorig2, icommcart
     $     ,mpiid,lflag,out)

      sor_del=delta
      ierr=isor_k
c Indirection is needed here because otherwise the finalize call
c seems to cause the return to fail. Probably unnecessary.
c      mpiid=myid

      end

c**********************************************************************
c Do a single relaxation.
      subroutine sorrelax(isor_k,Li,Lj,ni,nj,
     $     u,q,a,b,c,d,f,g,out,
     $     relax,rdelta,bcphi)

      integer isor_k,Li,Lj,ni,nj,bcphi
      logical out
      real u(Li,nj),q(Li,nj)
      real a(Li),b(Li),c(Li,nj),d(Li,nj),f(Li,nj),g(Lj,5)
      parameter (nd2=4)

      real relax,delta,rdelta

      rdelta=0.
c Hard-wired red-black relaxation order. isor_k odd => block is odd;
c odd meaning that (1,1) and hence (2,2) is active.
      km=mod(isor_k+1,2)
c (Remember that ni=rshield+1)
      do j=2,nj-1,1
         do i=2+mod(j+km,2),ni-1,2
        
c     Handle BC
            if(out.and.i.eq.ni-1) then
               delta=g(j,1)*u(i-1,j)+g(j,2)*u(i,j-1)+
     $              g(j,3)*u(i,j+1)+g(j,4)+g(j,5)*u(i,j)
               u(i,j)=u(i,j)+relax*delta
               delta=abs(delta)
                  if(delta.gt.rdelta)rdelta=delta
c     Handle the bulk calculations
            else
               expu=exp(u(i,j))
               dnum= a(i)*u(i+1,j)+b(i)*u(i-1,j) + c(i,j)
     $              *u(i,j+1)+d(i,j)*u(i,j-1) -f(i,j)*u(i,j)
     $              + q(i,j) - expu
               dden=f(i,j) + expu
               delta=relax*dnum/dden
               u(i,j)=u(i,j)+delta
               delta=abs(delta)
               if(delta.gt.rdelta)rdelta=delta
            endif
         enddo
      enddo

      end

c***********************************************************************
c The challenge here is to ensure that all processes decide to end
c at the same time. If not then a process will hang waiting for message.
c So we have to universalize the convergence test. All block must be
c converged, but the total spread depends on multiple blocks.
      subroutine testifconverged(sor_eps,delta,lconverged,sor_comm)
      include 'mpif.h'
      logical lconverged
      real delta,deltaR,sor_eps
      integer sor_comm,id
      
c Here we need to allreduce the data, selecting the maximum values,
c doing it in place.

c Loki's version of mpich seems not to like MPI_IN_PLACE
      call MPI_ALLREDUCE(delta,deltaR,1,MPI_REAL,MPI_MAX
     $     ,sor_comm,ierr)

      delta=deltaR

      if(delta.lt.sor_eps) then
         lconverged=.true.
      else
         lconverged=.false.
      endif

      end
c***********************************************************************
c Cut here and throw the rest away for the basic routines
c***********************************************************************
      subroutine fcalc_sor(Li,Lj,ni,nj,jac,mi,dconverge,k,sor_comm,myid2
     $     )


      include 'piccom.f'
      integer Li,Lj,ni,nj
      integer sor_comm,myid2
      integer nd
      parameter (nd=2,nd2=nd*2)

      real sor_jac,sor_eps,sor_del,jac,dconverge
      integer isor_k,isor_mi,mi
      

      common /sor2dctl/isor_mi,sor_jac,sor_eps,sor_del,isor_k

c testing arrays
      integer iuds(nd),ifull(nd)

      isor_mi=mi
      sor_jac=jac
      sor_eps=dconverge

      ifull(1)=Li
      ifull(2)=Lj
      iuds(1)=ni
      iuds(2)=nj
      
      ictl=1
  
      call sor2dmpi(sor_comm,Li,Lj,ni,nj,debyelen,bcphi, phi(1,0) ,rho(1
     $     ,0),ictl,ierr,myid2,idim1,idim2, apc(1),bpc(1),cpc(1,0)
     $     ,dpc(1,0),fpc(1,0),gpc(0,1))
      
      
      k=isor_k
      
      end
