c*********************************************************************
      function omlux(Ui,x)
c The functional dependence of OML collection vs drift U, and potl x.
      U=Ui
c More accurate to make the cut-off lower 1.e-2 vs 1.e-5
      if(U.lt.1.e-2) U=1.e-2
      omlux=( (1.+(1.+2.*x)/(2.*U**2))*(1.-erfcc(U))
     $     + exp(-U**2)/(sqrt(3.14159)*U) )*U
      end
c*****************************************************************
      function omlfloat(U,ZTei,Zmepi)
c Floating potential according to OML approximation.
c U is ratio of drift to ion thermal velocity
c ZTei is ZT_e/T_i
c Zmepi is sqrt(Z m_e 4\pi/m_i)
c Value returned is in units of ZT_e
      x1=2.*ZTei
 1    x=x1
      x1=-ZTei*alog( 0.25*(Zmepi/sqrt(ZTei))*omlux(U,x) )
c      write(*,*)x,omlux(u,x)
      if(abs(x1-x)/x.gt.1.e-5) goto 1
      omlfloat=x/ZTei
      end
c******************************************************************
      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
c*********************************************************************
      integer npts,ncur
      parameter (npts=50,ncur=3)
      real foml2d(npts,ncur),rmtoz(npts)
      real approx(npts,ncur)
      character*30 str1

      str1=' T!di!d='
      call pltinit(0.,float(npts),-5.,0.)
      call scalewn(1.,float(npts),-5.,0.,.true.,.false.)
      call axis()
      call axlabels('M','Floating Potential')
      emax=0.
      dmax=0.
      do j=1,ncur
         ZTei=10.**(j-1.)
         do i=1,npts
            rmtoz(i)=i
            Zmepi=sqrt(4.*3.1415926/(1837*rmtoz(i)))
            foml2d(i,j)=-omlfloat(0.,ZTei,Zmepi)
            x=1837.*rmtoz(i)/Ztei
c This approx is correct to better than 2% for Ti/Te <= 1.
            approx(i,j)=-0.36*log(x)
c     $           + .043*log10(Ztei/10.)*log(rmtoz(i)/80.)
     $           + .018*log(Ztei/10.)*log(rmtoz(i)/80.)
            write(*,*)rmtoz(i),foml2d(i,j),approx(i,j)
     $           ,foml2d(i,j)/approx(i,j)-1.
            e=abs(foml2d(i,j)/approx(i,j)-1.)
            if(e.gt.emax)emax=e
            d=abs(foml2d(i,j)-approx(i,j))
            if(d.gt.dmax)dmax=d
         enddo
      write(*,*)'Max frac error=',emax,'  Max difference=',dmax
      call color(j)
      call polyline(rmtoz,foml2d(1,j),npts)
      call dashset(1)
      call polyline(rmtoz,approx(1,j),npts)
      call dashset(0)
      call fwrite(1./ZTei,iwd,3,str1(9:))
      call legendline(.03,.02+.05*j,0,str1)
      enddo

      call pltend()
      end
