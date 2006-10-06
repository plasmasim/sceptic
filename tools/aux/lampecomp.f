      integer npts
      parameter (npts=20)
      real x(npts),y(npts),z(npts),z2(npts),y2(npts),z3(npts)
      real v(npts),f(npts)

      character*10 str1
c Data from Lampe et al POP03 for Ti=.01, a/lambda_s=.015
      parameter (nflampe=16)
      real flampe(nflampe),colflampe(nflampe),flm(nflampe),clm(nflampe)
      data flampe/0.,-25,-127.,-265.,-495.,-715.,-905.,-1195.,-1585.,
     $     -1925.,-2235.,-2535.,-2865.,-3215.,-3585.,-3955./
      data colflampe/0.,10.,50.,107.,207.,327.,437.,647.,967.,
     $     1327.,1697.,2177.,2787.,3567.,4537.,5637./
c upper point scaled positions and values:
      data cl1/5997./fl1/-3585./clv1/0.5/flv1/2./
c Dummy data
c      Ti=.01
      Ti=.1
      debyelen=666.7
      foml=9.67
      rmtoz=20.

      do kl=1,nflampe
c     Scale into Lampe's units:
         flampe(kl)=1.+flampe(kl)*flv1/fl1
         colflampe(kl)=colflampe(kl)*clv1/cl1
c         write(*,*)colflampe(kl),flampe(kl)
c     Scale into my units
         flm(kl)=flampe(kl)*foml
c     Unfortunately, Lampe uses lambda_s instead of lambda_De
         clm(kl)=colflampe(kl)*sqrt(2.*Ti)
     $        *sqrt(1+1./Ti)/debyelen
c         write(*,*)clm(kl),flm(kl)
      enddo
c The following form fits the original lampe results with astonishing 
c accuracy.
      do i=1,npts
         x(i)=0.5*(i-1)/(npts-1)
         y(i)=1.+log(1.+17.*x(i))
      enddo
      call autoplot(colflampe,flampe,nflampe)
      call dashset(1)
      call polyline(x,y,npts)
      call pltend()
c      Ti=.005
      do i=1,npts
         v(i)=10**(-4.+4.*(i-1)/(npts-1))
         f(i)=changlaf(Ti,-v(i),.1)
      enddo
      call fwrite(Ti,iwd,3,str1)
      call lautoplot(v,f,npts,.true.,.true.)
      call axlabels('Minus Probe Potential',
     $     'Chang Flux at Ti='//str1(1:iwd))
      call pltend()
      do i=1,npts
         x(i)=10**(-4.+6.*(i-1)/(npts-1))
         y(i)=changfloat(Ti,rmtoz,x(i),fluxi)
c         y2(i)= -3.4*log10(1.+sqrt(rmtoz/40.)*x(i)/.008)/log10(1./.008)
c compacter form:
c         y2(i)=-0.704*log(1.+19.76*sqrt(rmtoz)*x(i))
c Rounded:
c         y2(i)=-0.7*log(1.+20.*sqrt(rmtoz)*(x(i)))
c Good to higher collisionality:
         y2(i)=-0.52*(log(1.1+30.*sqrt(rmtoz)*x(i)))**1.15
c Crude low collisionality plateau:
c         y2(i)=min(y2(i),-(Ti)**.73*2.65)
c         z(i)=min(fluxi,100.)
         z(i)=fluxi
c This is a pretty good fit to the exact changlaframboise.
         z2(i)=-y2(i)/x(i)
c This version is even better:
         z3(i)=exp(y2(i))/sqrt(2.*3.1415926/(1837.*rmtoz))
c Reduced expression is the same only for old version:
c         z3(i)=(1.+19.76*sqrt(rmtoz)*x(i))**(-.704)
c     $        /sqrt(2.*3.1415926/(1837.*rmtoz))
      enddo
      call dashset(1)
      call lautoplot(x,y2,npts,.true.,.false.)
      call axlabels('collision frequency','floating potential')
      call legendline(.1,.15,0,' Fit')
      call dashset(0)
      call legendline(.1,.2,0,' Changfloat')
      call polyline(x,y,npts)
      call pltend()
      call lautoplot(x,z,npts,.true.,.false.)
      call axlabels('collision frequency','flux')
      call dashset(2)
      call polyline(x,z2,npts)
      call legendline(.1,.1,0,'z2')
      call dashset(3)
      call polyline(x,z3,npts)
      call legendline(.1,.15,0,'z3')
      call pltend()
      end
c********************************************************************
c Chang and Laframboise collisional ion flux:
c Needed to fix vti ambiguity there's a factor of order unity.
c First guess was low by factor 1.329:
      function changlaf(Ti,Vprobe,colnwt)
      kappai=sqrt(3.141592*2.*Ti)/colnwt
      expion=exp(Vprobe/Ti)
      fluxii=(1.+kappai)*(-Vprobe)/
     $     (Ti*(1.-(1.+kappai*Vprobe/Ti)*expion))
      changlaf=fluxii*Ti/(1+kappai)/(colnwt)
c net: flux=((-V)/colnwt) /(1-(1+k.V/Ti)*exp(V/Ti))
      end
c********************************************************************
      function changfloat(Ti,rmtoz,colnwt,fluxi)
      Vprobe=-2.
      c1=sqrt(2.*3.1415926/(rmtoz*1837.))
      do i=1,40
         fluxi=changlaf(Ti,Vprobe,colnwt)
c exp(Vprobe)\sqrt(Te/2\pi me)=fluxi \sqrt(Te/mi)
c Less convergent approaches.
c         Vnew=log(fluxi*c1)
c         dv=-(Vprobe-log(fluxi*c1))/(1.+5./(fluxi*c1))
         ev=exp(Vprobe)
         dv=-(ev-fluxi*c1)/(ev-fluxi*c1/Vprobe)
         Vnew=Vprobe+(.2+.8*(colnwt/(.01+colnwt)))*dv
         if(Vnew.ge.0.)Vnew=-1.e-20
         if(abs((Vnew-Vprobe)/Vprobe).lt.1.e-5)goto 100
         write(*,*)i,Vprobe,Vnew,fluxi
         Vprobe=Vnew
c         Vprobe=.5*Vnew+.5*Vprobe
c         if(.not.Vprobe.lt.0.)goto 110
      enddo
      write(*,*)'changfloat unconverged. Vprobe=',Vprobe
 100  continue
      changfloat=Vnew
      write(*,'(''changfloat: colnwt,V,it'',f8.5,f9.5,i3)')
     $     colnwt,Vprobe,i
      return
c Error case.
 110  continue
      changfloat=100
      return
      end
