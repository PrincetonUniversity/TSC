!     * subroutine axm2d2 *
!     *********************
!
!     purpose:
!        determine (xma,yma) of magnetic axis by polynomial fit
!
!        2d polynomial fit:   psi:= pfx(i)*xsi**i + pfy(j)*ynu**j +
!                                     +pfxy(i,j)*(xsi**i)*(ynu**j)
!
!                          using finite differences of second order
!
!     arguments:
!        psi(ndim,*)   ....2D scalar function                 IN
!        ii,jj         ....first guess of magn.axis indices   IN
!        xii,yjj       ....coord.of first guess: x(ii),y(jj)  IN
!        xma,yma       ....coordinates of magn.axis           OUT
!        pds (1)       ....psi(xma,yma)                       OUT
!        pds (2)       ....d/dx psi(xma,yma)                  OUT
!        pds (3)       ....d/dy psi(xma,yma)                  OUT
!        pds (4)       ....d2/dxdy psi(xma,yma)               OUT
!        pds (5)       ....d2/dxdx psi(xma,yma)               OUT
 
!        pds (6)       ....d2/dydy psi(xma,yma)               OUT
!        pfx (*)       ....coefficients                       OUT
!        pfy (*)       ....coefficients                       OUT
!        pfxy(ndima,*) ....coefficients                       OUT
!        dx,dy         ....meshsize in x- & y-direction       IN
!        isym          ....symmetry switch                    IN
!                          isym= 0 ...no symmetry
!                          isym= 1 ...symmetry around yjj
!
!     ------------------------------------------------------------------
!     SOURCE:    WOS:baldur11.axmFIT(axm2d2)                  25-08-86
!     LOAD:      wos:LIBWMAP        (axm2d2)   <cray>         25-08-86
!     ------------------------------------------------------------------
!                                                          wos   06/85
!***********************************************************************
!
      subroutine axm2d2 (psi,ndim,ii,jj,xii,yjj,x0,y0,xma,yma,pds,       &  
     &                       pfx,pfy,pfxy,ndima,dx,dy,isym)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ndim,ii,jj,ndima,isym,nout,ifct,jfct,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xii,yjj,x0,y0,xma,yma,pds,pfx,pfy,pfxy,dx,dy,psi,fsdx
      REAL*8 fsdxx,fsdy,fsdyy,pf00,xsima,ynuma,xsi,ynu
!============
      dimension  psi(ndim,*),pds(*)
      dimension  pfx(*),pfy(*),pfxy(ndima,*)
      data       nout/6/
!============      
!.990 format (1x,'axm2d2.cpu (sec):',3f12.6)
!     ------------------------------------------------------------------
      integer    S
!     ------------------------------------------------------------------
!
!ll   statement functions:
!l    --------------------
!
!l    second order:
         fsdx (ifct,jfct)= (      psi(ifct+1,jfct)                       &  
     &                       -    psi(ifct-1,jfct)) / (2.0_R8*dx)
         fsdxx(ifct,jfct)= (      psi(ifct+1,jfct)                       &  
     &                       -2.0_R8*psi(ifct  ,jfct)                    &  
     &                       +    psi(ifct-1,jfct)) / (dx*dx)
!
         fsdy (ifct,jfct)= (      psi(ifct,jfct+1)                       &  
     &                       -    psi(ifct,jfct-S)) / (2.0_R8*dy)
         fsdyy(ifct,jfct)= (      psi(ifct,jfct+1)                       &  
     &                       -2.0_R8*psi(ifct,jfct  )                    &  
     &                       +    psi(ifct,jfct-S)) / (dy*dy)
!     ------------------------------------------------------------------
!
!.....timing: ..........................
!...  call timeused (icpu0,io0,isys0)
!.......................................
!
!ll   calculate pf:(pfx,pfy,pfxy) directly:    <second  order>
!.....call printx ('(d) direct 2d  second  order  fit:$',2)
         S = 1
      if (isym.eq. 1)  S = -1
         pf00= psi(ii,jj)
         i = ii
         j = jj
         pfx(1)= fsdx(i,j)
         pfx(2)= fsdxx(i,j)/2.0_R8
         pfy(1)= fsdy(i,j)
         pfy(2)= fsdyy(i,j)/2.0_R8
!
          i = ii
          j = jj
         pfxy (1,1)      = (      fsdy (i +1,j )                         &  
     &                       -    fsdy (i -1,j )) / (2.0_R8*dx)
!.......................................................................
!
!ll   determine xma,yma:
         xsima= -pfx(1)/(2.0_R8*pfx(2))
         xma  = xii+xsima
         ynuma= -pfy(1)/(2.0_R8*pfy(2))
         yma  = yjj+ynuma
!.....timing: ..........................
!...  call timeused (icpu1,io1,isys1)
!...     cpu1= (icpu1-icpu0)*1.0e-06
!.......................................
!
!ll   compute pds; psima and derivatives:
!__      nor= 2
!__   call psi2df (xsima,ynuma,pf00,pfx,pfy,pfxy,ndima,nor,pds)
!.....direct:
      xsi = xsima
      ynu = ynuma
!
         pds(1)= pf00+pfx(1)*xsi+pfy(1)*ynu + pfxy(1,1)*xsi*ynu
         pds(2)=     +pfx(1)                + pfxy(1,1)    *ynu
         pds(3)=                +pfy(1)     + pfxy(1,1)*xsi
         pds(4)=                              pfxy(1,1)
         pds(5)= 0.0_R8
         pds(6)= 0.0_R8
!     ..................................................................
!
!.....timing: ..........................
!...  call timeused (icpu2,io2,isys2)
!...     cpu2= (icpu2-icpu1)*1.0e-06
!...     cpusum= cpu1+cpu2
!...  write (nout,990)  cpu1,cpu2,cpusum
!.......................................
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
