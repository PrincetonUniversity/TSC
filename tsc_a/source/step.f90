      subroutine step(x,y,neqn,h,eps,wt,start,                           &  
     &hold,k,kold,crash,phi,p,yp,psi)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER neqn,k,kold,l,ifail,kp1,kp2,km1,km2,ns,nsp1,i,im1,iq
      INTEGER nsm2,j,limit1,nsp2,limit2,ip1,knew
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,h,eps,wt,hold,phi,p,yp,psi,x,alpha,beta,sig,w,v,g
      REAL*8 gstr,two,twou,fouru,p5eps,round,absh,realns,temp1
      REAL*8 temp2,reali,temp3,temp4,temp5,temp6,tau,xold,erkm2
      REAL*8 erkm1,erk,rho,erkp1,hnew,r
      REAL*8 sum, err
!============
      logical start,crash,phase1,nornd
      dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
      dimension alpha(12),beta(12),sig(13),w(12),v(12),g(13),            &  
     &gstr(13),two(13)
      data twou,fouru/1.42E-15_R8,2.84E-15_R8/
      data two/2.0_R8,4.0_R8,8.0_R8,16.0_R8,32.0_R8,64.0_R8,128.0_R8,    &  
     & 256.0_R8,512.0_R8,                                                &  
     & 1024_R8,2048.0_R8,4096.0_R8,8192.0_R8/
      data gstr/.50_R8,.08330_R8,.04170_R8,.02640_R8,.01880_R8,          &  
     & .01430_R8,                                                        &  
     &.01140_R8,.009360_R8,.007890_R8,.006790_R8,.005920_R8,.005240_R8,  &  
     & .004680_R8/
      data g(1),g(2)/1.0_R8,.50_R8/,sig(1)/1.0_R8/
!============      
      crash=.true.
      if(abs(h).ge.fouru*abs(x)) goto 5
      h=sign(fouru*abs(x),h)
      return
5     p5eps=.50_R8*eps
      round=0.0_R8
      do 10 l=1,neqn
10    round=round+(y(l)/wt(l))**2
      round=twou*sqrt(round)
      if(p5eps.ge.round) goto 15
      eps=2.0_R8*round*(1.0_R8+fouru)
      return
15    crash=.false.
      if(.not.start) goto 99
      call der(x,y,yp,neqn)
      sum=0.0_R8
      do 20 l=1,neqn
      phi(l,1)=yp(l)
      phi(l,2)=0.0_R8
20    sum=sum+(yp(l)/wt(l))**2
      sum=sqrt(sum)
      absh=abs(h)
      if(eps.lt.16.0_R8*sum*h*h) absh=.250_R8*sqrt(eps/sum)
      h=sign(max(absh,fouru*abs(x)),h)
      hold=0.0_R8
      k=1
      kold=0
      start=.false.
      phase1=.true.
      nornd=.true.
      if(p5eps.gt.100.0_R8*round) goto 99
      nornd=.false.
      do 25 l=1,neqn
25    phi(l,15)=0.0_R8
99    ifail=0
100   kp1=k+1
      kp2=k+2
      km1=k-1
      km2=k-2
      if(h.ne.hold) ns=0
      ns=min(ns+1,kold+1)
      nsp1=ns+1
      if(k.lt.ns) goto 199
      beta(ns)=1.0_R8
      realns=ns
      alpha(ns)=1.0_R8/realns
      temp1=h*realns
      sig(nsp1)=1._R8
      if(k.lt.nsp1) goto 110
      do 105 i=nsp1,k
      im1=i-1
      temp2=psi(im1)
      psi(im1)=temp1
      beta(i)=beta(im1)*psi(im1)/temp2
      temp1=temp2+h
      alpha(i)=h/temp1
      reali=i
105   sig(i+1)=reali*alpha(i)*sig(i)
110   psi(k)=temp1
      if(ns.gt.1) goto 120
      do 115 iq=1,k
      temp3=iq*(iq+1)
      v(iq)=1.0_R8/temp3
115   w(iq)=v(iq)
      goto 140
120   if(k.le.kold) goto 130
      temp4=k*kp1
      v(k)=1.0_R8/temp4
      nsm2=ns-2
      if(nsm2.lt.1) goto 130
      do 125 j=1,nsm2
      i=k-j
125   v(i)=v(i)-alpha(j+1)*v(i+1)
130   limit1=kp1-ns
      temp5=alpha(ns)
      do 135 iq=1,limit1
      v(iq)=v(iq)-temp5*v(iq+1)
135   w(iq)=v(iq)
      g(nsp1)=w(1)
140   nsp2=ns+2
      if(kp1.lt.nsp2) goto 199
      do 150 i=nsp2,kp1
      limit2=kp2-i
      temp6=alpha(i-1)
      do 145 iq=1,limit2
145   w(iq)=w(iq)-temp6*w(iq+1)
150   g(i)=w(1)
199   continue
      if(k.lt.nsp1) goto 215
      do 210 i=nsp1,k
      temp1=beta(i)
      do 205 l=1,neqn
205   phi(l,i)=temp1*phi(l,i)
210   continue
215   do 220 l=1,neqn
      phi(l,kp2)=phi(l,kp1)
      phi(l,kp1)=0.0_R8
220   p(l)=0.0_R8
      do 230 j=1,k
      i=kp1-j
      ip1=i+1
      temp2=g(i)
      do 225 l=1,neqn
      p(l)=p(l)+temp2*phi(l,i)
225   phi(l,i)=phi(l,i)+phi(l,ip1)
230   continue
      if(nornd) goto 240
      do 235 l=1,neqn
      tau=h*p(l)-phi(l,15)
      p(l)=y(l)+tau
235   phi(l,16)=(p(l)-y(l))-tau
      goto 250
240   do 245 l=1,neqn
245   p(l)=y(l)+h*p(l)
250   xold=x
      x=x+h
      absh=abs(h)
      call der(x,p,yp,neqn)
      erkm2=0.0_R8
      erkm1=0.0_R8
      erk=0.0_R8
      do 265 l=1,neqn
      temp3=1.0_R8/wt(l)
      temp4=yp(l)-phi(l,1)
      if(km2) 265,260,255
255   erkm2=erkm2+((phi(l,km1)+temp4)*temp3)**2
260   erkm1=erkm1+((phi(l,k)+temp4)*temp3)**2
265   erk=erk+(temp4*temp3)**2
      if(km2) 280,275,270
270   erkm2=absh*sig(km1)*gstr(km2)*sqrt(erkm2)
275   erkm1=absh*sig(k)*gstr(km1)*sqrt(erkm1)
280   temp5=absh*sqrt(erk)
      err=temp5*(g(k)-g(kp1))
      erk=temp5*sig(kp1)*gstr(k)
      knew=k
      if(km2) 299,290,285
285   if(max(erkm1,erkm2).le.erk) knew=km1
      goto 299
290   if(erkm1.le..50_R8*erk) knew=km1
299   if(err.le.eps) goto 400
      phase1=.false.
      x=xold
      do 310 i=1,k
      temp1=1.0_R8/beta(i)
      ip1=i+1
      do 305 l=1,neqn
305   phi(l,i)=temp1*(phi(l,i)-phi(l,ip1))
310   continue
      if(k.lt.2) goto 320
      do 315 i=2,k
315   psi(i-1)=psi(i)-h
320   ifail=ifail+1
      temp2=.5_R8
      if(ifail-3) 335,330,325
325   if(p5eps.lt..250_R8*erk) temp2=sqrt(p5eps/erk)
330   knew=1
335   h=temp2*h
      k=knew
      if(abs(h).ge.fouru*abs(x)) goto 340
      crash=.true.
      h=sign(fouru*abs(x),h)
      eps=eps+eps
      return
340   goto 100
400   kold=k
      hold=h
      temp1=h*g(kp1)
      if(nornd) goto 410
      do 405 l=1,neqn
      rho=temp1*(yp(l)-phi(l,1))-phi(l,16)
      y(l)=p(l)+rho
405   phi(l,15)=(y(l)-p(l))-rho
      goto 420
410   do 415 l=1,neqn
415   y(l)=p(l)+temp1*(yp(l)-phi(l,1))
420   call der(x,y,yp,neqn)
      do 425 l=1,neqn
      phi(l,kp1)=yp(l)-phi(l,1)
425   phi(l,kp2)=phi(l,kp1)-phi(l,kp2)
      do 435 i=1,k
      do 430 l=1,neqn
430   phi(l,i)=phi(l,i)+phi(l,kp1)
435   continue
      erkp1=0.0_R8
      if(knew.eq.km1.or.k.eq.12) phase1=.false.
      if(phase1) goto 450
      if(knew.eq.km1) goto 455
      if(kp1.gt.ns) goto 460
      do 440 l=1,neqn
440   erkp1=erkp1+(phi(l,kp2)/wt(l))**2
      erkp1=absh*gstr(kp1)*sqrt(erkp1)
      if(k.gt.1) goto 445
      if(erkp1.ge..50_R8*erk) goto 460
      goto 450
445   if(erkm1.le.min(erk,erkp1)) goto 455
      if(erkp1.ge.erk.or.k.eq.12) goto 460
450   k=kp1
      erk=erkp1
      goto 460
455   k=km1
      erk=erkm1
460   hnew=h+h
      if(phase1) goto 465
      if(p5eps.ge.erk*two(k+1)) goto 465
      hnew=h
      if(p5eps.ge.erk) goto 465
      temp2=k+1
      r=(p5eps/erk)**(1.0_R8/temp2)
      hnew=absh*max(.50_R8,min(.90_R8,r))
      hnew=sign(max(hnew,fouru*abs(x)),h)
465   h=hnew
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
