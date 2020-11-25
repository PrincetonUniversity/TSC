subroutine w20wguess ( wz )
  use w20data
  implicit none

  COMPLEX*16, intent(inout) :: wz

  COMPLEX*16 :: wz1, wz2, H_q, E1, E2
  COMPLEX*16 :: iu
  REAL*8     :: G_ave

  G_ave = 1.0
  IU = (0.0, 1.0)

  H_q = cmplx(0.0,0.5*shat/q)

  E1 = ftrt*(1.0+zflh)-0.5*gne+0.5*zflh*tauh*(gne+gti) &
     + (gm+ftrt)*G_ave+H_q*(1.0+ftr*tauh)

  E1 = E1 * 0.5 / ( 1.0 + zflh )


  E2 =  (  ( 0.5 * tauh * gm * ( gti - tvr*gne ) + bta )*( G_ave + H_q ) &
       - 0.5 * tauh * ftr * ( gni - zflh * tauh * (gti + gni))  ) &
       / ( 1.0 + zflh )

  wz1=-E1+SQRT(E1*E1-E2)
  wz2=-E1-SQRT(E1*E1-E2)

  wz=wz1

  IF(DIMAG(wz2).GT.DIMAG(wz1)) wz=wz2

  if ( dimag(wz) .le. 0.01 ) then
     wz = wz - iu*dimag(wz) + iu*0.01
  end if

end subroutine w20wguess
