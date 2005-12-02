      integer npartmax,npart,nr,nth,ndim,np
      parameter (npartmax=200000,nr=401,nth=101,np=1,ndim=6)
c NGP definitions
      logical LCIC
      integer NRUSED,NTHUSED,NRFULL,NTHFULL
      parameter (LCIC=.false.)
      parameter (NRUSED=nr-1,NTHUSED=nth-1,NRFULL=nr,NTHFULL=nth)
c Positions and velocities of particles (6-d phase-space).
      real xp(ndim,npartmax)
c The particle number
      real psum(0:NRFULL,0:NTHFULL)
c The sum of particle radial velocities
      real vrsum(0:NRFULL,0:NTHFULL)
c Sum of theta velocities
      real vtsum(0:NRFULL,0:NTHFULL)
c Sum of phi velocities
      real vpsum(0:NRFULL,0:NTHFULL)
c The sum of particle velocities squared
      real v2sum(0:NRFULL,0:NTHFULL)
c The sum of radial particle velocities squared
      real vr2sum(0:NRFULL,0:NTHFULL)
c The sum of non-radial particle velocities squared
      real vtp2sum(0:NRFULL,0:NTHFULL)
c The potential normalized to Te/e
      real phi(0:NRFULL,0:NTHFULL)
c Charge density
      real rho(0:NRFULL,0:NTHFULL)
c Mag field
      real Bz

      real pi
      parameter (pi=3.1415927)
      real cerr,bdyfc,Ti,vd
      logical diags,lplot
      integer myid,numprocs
      common /piccom/xp,npart,
     $     psum,vrsum,vtsum,vpsum,v2sum,vr2sum,vtp2sum,
     $     phi,rho,cerr,bdyfc,Ti,vd,diags,lplot
     $     ,myid,numprocs,Bz

c Radius mesh
      real r(0:NRFULL),rcc(0:NRFULL)
c Theta angle mesh
      real th(0:NTHFULL),tcc(0:NTHFULL)
c Theta mesh radians
      real thang(0:NTHFULL)
c Inverse of volume of radial shells.
      real volinv(0:NRFULL)
c Precalculation functions
      integer nrpre,ntpre
      parameter (nrpre=1000,ntpre=100)
      integer irpre(nrpre),itpre(ntpre)
      real rfac,tfac
c Non-uniform handling quantities.
      real hr(0:nr+1),zeta(0:nr+1),zetahalf(0:nr+1),
     $     cminus(nr),cmid(nr),cplus(nr)
c Lower limit of averaging range. 0.6 by default
      real avelim

      common /meshcom/r,rcc,th,tcc,thang,volinv,irpre,itpre,rfac,tfac,
     $     hr,zeta,zetahalf,cminus,cmid,cplus,avelim

      integer nvel
      parameter (nvel=50)
      real Qcom(nth)
      real Gcom(nvel,nth)
      real Vcom(nvel)
      common /rancom/Gcom,Vcom,Qcom

      integer nvmax,nrein,ninner,nstepmax
      parameter (nvmax=60,nstepmax=4000)
      real nvdiag(nvmax),nvdiagave(nvmax),vdiag(nvmax)
      real vrange
      real diagrho(nr),diagphi(nr)
      real diagvr(nr,nth)
c Total particle flux to probe
      real fluxprobe(nstepmax)
c Number of particles striking probe in theta cell
      integer ninthstep(nth,nstepmax)
      integer ninth(nth)
c Number of particles reinjected per theta cell.
      integer noutrein(nth),ivoutrein(nth)
      real rhoinf
      common /diagcom/nvdiag,nvdiagave,vdiag,vrange,diagrho,diagphi,
     $     nrein,ninner,fluxprobe,ninthstep,ninth,rhoinf,noutrein
     $     ,ivoutrein,diagvr



