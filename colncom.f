c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c This code is copyright(c) (2003-5) Ian H Hutchinson hutch@psfc.mit.edu
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build
c or run the code, you do so at your own risk.
c
c Version 2.6 Aug 2005.
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c Collision data common.
c Drift velocity of neutral distribution
      real vneutral
c Temperature of neutral distribution
      real Tneutral
c Required acceleration to give rise to specified ion average drift.
      real Eneutral
c Number of collisions this time step
      integer NCneutral
      common /colcom/vneutral,Tneutral,Eneutral,NCneutral
