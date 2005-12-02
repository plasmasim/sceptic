c___________________________________________________________________________
c
c This code is copyright (c) Ian H Hutchinson (2004) hutch@psfc.mit.edu.
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the paper I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002).  The code may not be
c redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build or run
c the code, you do so at your own risk.
c___________________________________________________________________________
c**************************************************************************
c Calculate the normalized flux density to a Coulomb sphere collecting ions
c in a shifted Maxwellian distribution.
c
c I.H.Hutchinson Sep 2002.
c Version using tricks with Bessel functions to avoid overflows 13 Sep 02.

c Return is flux density normalized to the stationary maxwellian one-way flux.
      real function coulflux(betapass,chipass,vdpass)
c Angle between the surface element normal and the drift direction, radians.
      real betapass
c Minus Potential on the sphere normalized to the Maxwellian Temperature.
      real chipass
c Drift velocity normalized to sqrt(T/m).
      real vdpass
      external uintegrand
      real umax,umin
      common /coulparam/cosbeta,sinbeta,chi,vd

c Convert to normalization by sqrt(2T/m) for calculation.
      vd=vdpass/sqrt(2.)
      umin=1.e-7
      umax=max(5.,abs(vd)+4.)
c Formulas are written for angle pi-beta, i.e. the inward normal.
c So convert to that.
      cosbeta=-cos(betapass)
      sinbeta=sin(betapass)
      chi=chipass
      call qrombu(uintegrand,umin,umax,coulflux)
      coulflux=coulflux*4.
      end
c*********************************************************************
      real function uintegrand(u)
      common /coulparam/cosbeta,sinbeta,chi,vd
      common /uiny/uval,tuu1,tuu2,ymax
      external sintg

      uval=u
      tuu1=2.*uval*vd*cosbeta
      tuu2=2.*uval*vd*sinbeta
      u2=u**2
      if(u2.le.chi*1.e-12)then
         ymax=1.e6
      else
         ymax=sqrt(1.+chi/u2)
      endif
c      ymin=1.e-8*ymax
      smin=1.e-10
      call qromby(sintg,smin,1.,sresult)
c      uintegrand=sresult*u2*u*exp(-(u2+vd**2))
c Take the exponential inside the inner integral to help prevent overflow.
      uintegrand=sresult*u2*u
      end
c*********************************************************************
      real function sintg(s)
c integrand for s-integration. y/ymax=s^2(3-2s).
      common /coulparam/cosbeta,sinbeta,chi,vd
      common /uiny/uval,tuu1,tuu2,ymax
      y=ymax*s**2*(3.-2.*s)
      ea1=tuu1*cosalpha(chi,y,uval)
      ea2=tuu2*sinalpha(chi,y,uval)
      sintg=exp(ea1+abs(ea2)-(uval**2+vd**2))*
     $     ebessi0(ea2)
     $     *y
     $     *ymax*6.*(1.-s)*s
      if(.not.sintg.le.1.e30)then
         write(*,*) 'Coulflux Overflow in integrand.'
         write(*,*) cosalpha(chi,y,uval),sinalpha(chi,y,uval),tuu1,tuu2
         stop
      endif
      end
c*********************************************************************
c*********************************************************************
      real function cosalpha(chi,y,u)
      real chi,y,u

      u2y=u**2*y
      if(u2y.le.chi*1.e-10)then
         cosalpha=1.
      else
         cbu=-chi/(2.*u2y)
         sarg=max(0.,1-2.*cbu*y-y**2)
         cosalpha=(sqrt(sarg)+(y+cbu)*cbu)/(1+cbu**2)
      endif
      end
c*********************************************************************
      real function sinalpha(chi,y,u)
      real chi,y,u

      u2y=u**2*y
      if(u2y.le.chi*1.e-10)then
         sinalpha=0.
      else
         cbu=-chi/(2.*u2y)
         sarg=max(0.,1-2.*cbu*y-y**2)
         sinalpha=(-sqrt(sarg)*cbu+(y+cbu))/(1+cbu**2)
      endif
      end
c**********************************************************************
c**********************************************************************
c Not used here.
      FUNCTION BESSI0(X)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
c**********************************************************************
c**********************************************************************
      FUNCTION eBESSI0(X)
c exp(-|x|)I_o(x)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      ax=abs(x)
      IF (Ax.LT.3.75) THEN
        Y=(X/3.75)**2
        eBESSI0=exp(-ax)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=3.75/AX
        eBESSI0=(1./SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
c**********************************************************************
c Based on Numerical Recipes.
      SUBROUTINE QROMBu(FUNC,A,B,SS)
      save
      external func
      PARAMETER (EPS=1.E-6, JMAX=12, JMAXP=JMAX+1, K=5, KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)

      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZDu(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINTu(H(J-KM),S(J-KM),K,0.,SS,DSS)
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
      write(*,*) 'Too many steps u.'
      END

      SUBROUTINE QROMBy(FUNC,A,B,SS)
      save
      external func
c jmax=12 is enough to get six decimal places of agreement.
      PARAMETER (EPS=1.E-5, JMAX=12, JMAXP=JMAX+1, K=5, KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)

      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZDy(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINTy(H(J-KM),S(J-KM),K,0.,SS,DSS)
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
      write(*,*) 'Too many steps y. Delta=',dss
      END

      SUBROUTINE POLINTu(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      SUBROUTINE POLINTy(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c**********************************************************************
      SUBROUTINE TRAPZDU(FUNC,A,B,S,N)
      save
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
c**********************************************************************
      SUBROUTINE TRAPZDY(FUNC,A,B,S,N)
      save
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

