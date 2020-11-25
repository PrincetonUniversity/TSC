module trdatbuf_iface

  ! interface for "public" routines in trdatbuf_lib
  ! also define constants useful for this interface...

  include 'trdatbuf_constants.incl'

  interface

     SUBROUTINE TDB_DATALO(d,ISIZE,ILOC)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d          ! data buffer
       integer, intent(in) :: isize  ! size of chunk to allocate
       integer, intent(out) :: iloc  ! address of chunk allocated
     end subroutine TDB_DATALO

     SUBROUTINE TDB_WORKALO(d,ISIZE,ILOC)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d          ! data buffer
       integer, intent(in) :: isize  ! size of chunk to allocate
       integer, intent(out) :: iloc  ! address of chunk allocated
     end subroutine TDB_WORKALO

     subroutine tdb_find_nmoms(d,imoms,icirc)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d          ! data buffer
       integer, intent(out) :: imoms,icirc
       ! return imoms = # of moments (0:imoms) in Fourier expansion
       ! return icirc = 1 if circular flux surfaces (minor & major radius only)
       !                are in use (this is very rare).
     end subroutine tdb_find_nmoms

     SUBROUTINE TDB_LOOKUP1(d,ztime,it1,it2,zf1,zf2)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       real*8, intent(in) :: ztime
       integer, intent(out) :: it1,it2  ! time offsets btw which ZTIME lies
       real*8, intent(out) :: zf1,zf2   ! linear interpolation factors
     end subroutine TDB_LOOKUP1

     real*8 function tdb_fintrp1(d,iloc,it1,it2,zf1,zf2)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       integer, intent(in) :: iloc     ! data location
       integer, intent(in) :: it1,it2  ! time offsets btw which ZTIME lies
       real*8, intent(in) :: zf1,zf2   ! linear interpolation factors
     end function TDB_FINTRP1

     subroutine tdb_rmp_bdy(d,ztime,zr1,zr2)
       !  find midplane boundary intercept radii at the indicated time

       use trdatbuf_obj_tsc
       implicit NONE

       type(trdatbuf) :: d          ! data buffer object
       real*8, intent(in) :: ztime  ! time (seconds)

       real*8, intent(out) :: zr1   ! inner intercept major radius (cm)
       real*8, intent(out) :: zr2   ! outer intercept major radius (cm)
     end subroutine tdb_rmp_bdy

     logical function tdb_logchk_special(d,zitem,iwarn)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       character*(*), intent(in) :: zitem   ! item queried
       integer, intent(out) :: iwarn        ! 0: OK; 1: item not recognized
     end function TDB_LOGCHK_SPECIAL

     logical function tdb_logchk_nbi(d,zsubset,iwarn)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       character*(*), intent(in) :: zsubset ! NBI subset data item queried
       integer, intent(out) :: iwarn        ! 0: OK; 1: item not recognized
     end function TDB_LOGCHK_NBI

     integer function tdb_nchan_find(d,z2char)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       character*(*), intent(in) :: z2char  ! "NB" or "EC" or "LH" or (IC)"RF"
     end function tdb_nchan_find

     subroutine tdb_find_next_sawtime(d,zprev_time,ztime)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       real*8, intent(in) :: zprev_time   ! preceding time (input)
       real*8, intent(out) :: ztime       ! time found or a huge time (output)
     end subroutine TDB_FIND_NEXT_SAWTIME

     subroutine tdb_num_sawtimes(d,ntsaw)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       integer, intent(out) :: ntsaw  ! number of sawtooth times 
     end subroutine TDB_NUM_SAWTIMES

     subroutine tdb_sawtimes(d,ztimes,inum)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       integer, intent(in) :: inum
       real*8, intent(out) :: ztimes(inum)
     end subroutine TDB_SAWTIMES

     subroutine tdb_post_sawtime(d,ztime_pre,ztime_post,iwarn)
       use trdatbuf_obj_tsc
       implicit NONE
       !  given event start (SAWTOOTH) time return event end time
       type (trdatbuf) :: d
       real*8, intent(in) :: ztime_pre   ! time at start of event
       real*8, intent(out) :: ztime_post ! time at end of event
       integer, intent(out) :: iwarn     ! 0=OK; 1=no sawtooth events
       ! 2= sawtooth events exist but ztime_pre
       !    does not match any of them...
     end subroutine TDB_POST_SAWTIME

     subroutine tdb_post_peltime(d,ztime_pre,ztime_post,iwarn)
       use trdatbuf_obj_tsc
       implicit NONE
       !  given event start (PELLET) time return event end time
       type (trdatbuf) :: d
       real*8, intent(in) :: ztime_pre   ! time at start of event
       real*8, intent(out) :: ztime_post ! time at end of event
       integer, intent(out) :: iwarn     ! 0=OK; 
                                    ! 2= no such event time found...
     end subroutine TDB_POST_PELTIME

     subroutine tdb_pwrdata_avg(d,znbi,ierr)
       use trdatbuf_obj_tsc
       use trdatbuf_aux
       implicit NONE
       !  get beam parameters averaged over a time step
       type (trdatbuf) :: d   ! trdat data object
       type (pwrget) :: znbi  ! specification of desired parameters; return val
       integer, intent(out) :: ierr   ! completion code (0=OK)
     end subroutine tdb_pwrdata_avg

     ! data on MHD equilibrium boundary...

     subroutine tdb_bdy_midplane(d,ztime,z_midp,rmaj_midp,rmin_midp)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       real*8, intent(in) :: ztime    ! time (seconds)
       real*8, intent(out) :: Z_midp  ! midplane offset from VV midplane (cm)
       real*8, intent(out) :: Rmaj_midp  ! midplane - plasma major radius (cm)
       real*8, intent(out) :: Rmin_midp  ! midplane - plasma minor radius or
       !                                   half-width (cm)
     end subroutine tdb_bdy_midplane

     subroutine tdb_bdy_timefac(d,ztime,iti,zfi)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       real*8, intent(in) :: ztime    ! time (seconds)
       integer, intent(out) :: iti    ! time bin index
       real*8, intent(out) :: zfi     ! time interpolation factor
       !                                   half-width (cm)
     end subroutine tdb_bdy_timefac

     real*8 function tdb_bdy_getmom(d,iti,zfi,imom,itype)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       integer, intent(in) :: iti    ! time bin index
       real*8, intent(in) :: zfi     ! time interpolation factor
       integer, intent(in) :: imom   ! moment index
       integer, intent(in) :: itype  ! moment type
       ! TBD_MOMS_RCOS or TBD_MOMS_RSIN or TBD_MOMS_ZCOS or TBD_MOMS_ZSIN
       ! ...function value in cm
     end function tdb_bdy_getmom
     
     ! data on impurities...

     logical function tdb_tokamakium(d)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d   ! return TRUE if "tokamakium" impurity is present
     end function TDB_TOKAMAKIUM

     subroutine tdb_imp_elems(d,xzimps,aimps)
       use trdatbuf_obj_tsc
       implicit NONE
       !  for multi-impurity model only...
       type (trdatbuf) :: d
       real*8, dimension(:) :: xzimps,aimps  ! Z & A vectors
     end subroutine tdb_imp_elems

     subroutine tdb_imp_ions(d,xzimpx,xzimpxs,aimpx)
       use trdatbuf_obj_tsc
       implicit NONE
       !  for multi-impurity model only...
       type (trdatbuf) :: d
       real*8, dimension(:) :: xzimpx !  charge on ion e.g. Z(C+4)=4
       real*8, dimension(:) :: xzimpxs !  atomic number e.g. Z(C+4)=6
       real*8, dimension(:) :: aimpx  !  atomic weight e.g. A(C+4)=14 if Carbon-14
     end subroutine tdb_imp_ions

     ! time invariant information in PH.DAT file...
     subroutine tdb_ccwchk(d,i_bccw,i_jccw)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       !  ccw means "counter-clockwise looking down from above"
       integer, intent(out) :: i_bccw ! 1: ccw, -1: cw, 0: unspecified (B_phi)
       integer, intent(out) :: i_jccw ! 1: ccw, -1: cw, 0: unspecified (J_phi)
     end subroutine tdb_ccwchk

     subroutine tdb_ripple0(d,icoils,zphi)
       use trdatbuf_obj_tsc
       implicit NONE
       type (trdatbuf) :: d
       !  information on TF ripple field
       integer, intent(out) :: icoils  ! no. of TF coils (TF ripples)
       real*8, intent(out) :: zphi ! shift of 1st coil rel. phi=0 (radians)
     end subroutine tdb_ripple0

     subroutine tdb_ki2mod(d,ifki2,imod)
       use trdatbuf_obj_tsc
       implicit NONE
       !
       !  return information on chi(i) model
       !
       type (trdatbuf) :: d
       integer, intent(out) :: ifki2 ! .gt.0 if chi(i) profile data is present
       integer, intent(out) :: imod  ! .gt.0 if data is multiplier on Neoclassical
       ! model; if .eq.0, the data is chi(i) itself, cm**2/sec...
     end subroutine tdb_ki2mod

     subroutine tdb_tlims(d,zd_tinit,zd_ftime)
       use trdatbuf_obj_tsc
       implicit NONE
       !
       !  return time range covered by trdatbuf data
       !
       type (trdatbuf) :: d
       real*8, intent(out) :: zd_tinit  ! start time
       real*8, intent(out) :: zd_ftime  ! stop time
     end subroutine tdb_tlims

     subroutine tdb_chkprof(d,istop)
       use trdatbuf_obj_tsc
       implicit NONE
       !
       !  check all input profile controls & x axes
       !  messages may be generated; all errors should result in cancellation
       !  of job.
       !
       type (trdatbuf) :: d
       integer, intent(out) :: istop    ! =0: OK; =1: error reported.
     end subroutine tdb_chkprof

      subroutine tdb_symini(d,nzones)
        !
        !  set up workspace for presymmetrization & mapping of profiles
        !  size of workspace depends on target application's number of zones
        !  and the number of zones in the input data
        !
        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        integer, intent(in) :: nzones   ! no. of radial zones in caller's grid
      end subroutine tdb_symini

      subroutine tdb_ntimes(d,intime1,intime2)
        !  return sizes of tdb input 1d f(t) timebase and 2d f(t,x) timebase
        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        integer, intent(out) :: intime1
        integer, intent(out) :: intime2
      end subroutine tdb_ntimes
      
      subroutine tdb_onoff_atime(d,z2char,zthresh,zdtfix,ton,toff,ierr)
        ! for power-channel data (NB,EC,LH,RF) define the on and off times:
        !   ton = min(on times for any channel) (*output*)
        !   toff = max(off times for any channel) (*output*)
        !
        !   z2char => chooses heating channel:
        !     "nb" or "NB" -- neutral beams
        !     "ec" -- ECH/ECCD
        !     "lh" -- Lower Hybrid
        !     "rf" -- ICRF   
        !   all tests of z2char value are case insensitive
        ! 
        !   zthresh = threshhold:
        !     if positive -- a power, in watts, must be < 20% of the maximum
        !                    power ever occurring on any channel
        !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
        !                    power threshold becomes -zthresh * (maximum power
        !                    ever occurring on any channel).
        !
        !   zdtfix -- time to search from first/last powers satisfying thresh-
        !             hold, for an actual 0 or negative power...
        !
        !   ierr is set only if there is no data or if z2char is unrecognized;
        !   if the zthresh limit has to be adjusted to conform to rules, a 
        !   warning message is written but ierr is not set.
        !
        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
        real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
        real*8, intent(in) :: zdtfix         ! time to search from threshhold
        real*8, intent(out) :: ton,toff      ! on/off times (seconds)
        integer, intent(out) :: ierr         ! completion code; 0=OK
      end subroutine tdb_onoff_atime

      subroutine tdb_onoff_times(d,z2char,zthresh,zdtfix,tonarr,toffarr,ierr)
        ! for power-channel data (NB,EC,LH,RF) define the on and off times:
        !   tonarr(i) = on time for channel #i (i.e. beam or antenna) (*output*)
        !   toffarr(i) = off time for channel #i (*output*)
        !
        !   z2char => chooses heating channel:
        !     "nb" or "NB" -- neutral beams
        !     "ec" -- ECH/ECCD
        !     "lh" -- Lower Hybrid
        !     "rf" -- ICRF   
        !   all tests of z2char value are case insensitive
        ! 
        !   zthresh = threshhold:
        !     if positive -- a power, in watts, must be < 20% of the maximum
        !                    power ever occurring on any channel
        !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
        !                    power threshold becomes -zthresh * (maximum power
        !                    ever occurring on any channel).
        !
        !   zdtfix -- time to search from first/last powers satisfying thresh-
        !             hold, for an actual 0 or negative power...
        !
        !   ierr is set only if there is no data or if z2char is unrecognized;
        !   if the zthresh limit has to be adjusted to conform to rules, a 
        !   warning message is written but ierr is not set.
        !
        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
        real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
        real*8, intent(in) :: zdtfix         ! time to search from threshhold
        real*8, intent(out), dimension(:) :: tonarr,toffarr  ! on/off times (seconds)
        integer, intent(out) :: ierr         ! completion code; 0=OK
      end subroutine tdb_onoff_times

      logical function tdb_fallback(d,ztri,ztime,zresult)
        !
        ! check for "fall back to time series" data
        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        character*(*), intent(in) :: ztri  ! trigraph of candidate data
        ! for "fall back to time series"
        real*8, intent(in) :: ztime        ! time to which to interpolate
        real*8, intent(out) :: zresult     ! result, if data is found
        ! the function value is .TRUE., and the interpolated value is returned,
        ! if "fallback to time series" data is present; otherwise the function
        ! returns .FALSE. and zresult=ZERO.
      end function tdb_fallback

      REAL*8 function tdb_get_rpld(d,zRi,zYi,ier)
        !
        ! get TF ripple magnitude in form: log(B~/B)
        !
        use trdatbuf_obj_tsc
        IMPLICIT NONE
        type (trdatbuf) :: d
        REAL*8,intent(in) :: zRi,zYi ! (R,Y) location, cm
        integer, intent(out) :: ier ! exit code, 0=OK
      end function tdb_get_rpld
 
      subroutine tdb_getmmx(d,zt,xib,zrmc2,zymc2,mj,mimom,lcentr,nzp1,zdrshaf)
        use trdatbuf_obj_tsc
        IMPLICIT NONE
 
        type (trdatbuf) :: d
        REAL*8,intent(in) :: zt    ! time at which to fetch moments
        integer, intent(in) :: mj  ! flux surface array dim. for zrmc2,zymc2
        integer, intent(in) :: mimom ! moments array dimension for zrmc2,zymc2

        real*8,intent(in) :: xib(mj) ! flux surfaces (xib(j)=sqrt(Phi/Philim) 
        ! at surface "j"); xib(lcentr:lcentr+nzp1-1) are used
  
        REAL*8,intent(out) :: zrmc2(mj,0:mimom,2),zymc2(mj,0:mimom,2)
        ! asymmetric Fourier moments set at xi bdys (lcentr:lcentr+nzp1-1)
 
        REAL*8,intent(out) :: zdrshaf(mj)
        ! "Shafranov" shift of interior surfaces
        ! relative to the boundary surface
 
        integer :: lcentr        ! index to magnetic axis
        integer :: nzp1          ! no. of surfaces including mag. axis
      end subroutine tdb_getmmx

      real*8 function tdb_zyget(d,zid,i)

        !  multiple impurity model (sim) helper function -- get charge state
        !  or get non-spatial / non-temporal 3rd coordinate for general 3d
        !  profile input...

        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        character*3,intent(in) :: zid  ! profile set identifier e.g. "SIM"
        integer,intent(in) :: i        ! get i'th element

      end function tdb_zyget
     
      subroutine tdb_profin(d,t,ierr)

        ! retrieve profile from trdatbuf buffer (d)
        ! to target time & spatical grid as specified in (t),
        ! desired profile also identified in (t)
        ! return ierr=0 if successful.

        use trdatbuf_obj_tsc
        use trdatbuf_aux
        implicit NONE
        type (trdatbuf) :: d
        type (profget) :: t
        integer, intent(out) :: ierr

      end subroutine tdb_profin

      subroutine tdb_maksym(d,ierr)

        ! presymmetrize all pre-symmetrizable profile data
        ! over the time range available

        use trdatbuf_obj_tsc
        implicit NONE
        type (trdatbuf) :: d
        integer, intent(out) :: ierr

      end subroutine tdb_maksym


     !! interface trdatbuf_write
     !! end interface trdatbuf_write

     !! interface trdatbuf_read_tsc
     !! end interface trdatbuf_read_tsc

     subroutine tdb_nzones(nzones,nrmaj,nrsym)
       integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)
       integer, intent(out) :: nrmaj  ! (2*nzones+1) no. of flux surface
       !                                midplane intercepts
       integer, intent(out) :: nrsym  ! (4*nzones+5) size of midplane test
       !                                grid (zone bdys & surfaces + 2
       !                                extrapolation points at each end).
     end subroutine tdb_nzones

     subroutine tdb_xilmp(nzones,xibdys,xilmp)
       integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)
       real*8, dimension(:), intent(in) :: xibdys  ! (nzones+1) x @ zone bdys
       real*8, dimension(:), intent(out) :: xilmp  ! (2*nzones+1) x @midplane
       !                                flux surface intercepts
     end subroutine tdb_xilmp

     subroutine tdb_xisymp(nzones,xibdys,rmajmp,xirsym,rmjsym)
       integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)
       real*8, dimension(:), intent(in) :: xibdys  ! (nzones+1) x @ zone bdys
       real*8, dimension(:), intent(in) :: rmajmp  ! (2*nzones+1) major radius
       !         at flux surface midplane intercepts (cm)
       real*8, dimension(:), intent(out) :: xirsym ! (4*nzones+5) x @midplane
       !         (double-density grid with 2 pt extension @ bdys)
       real*8, dimension(:), intent(out) :: rmjsym ! (4*nzones+5) R @midplane
       !         (double-density grid with 2 pt extension @ bdys) (cm)
     end subroutine tdb_xisymp

      subroutine tdb_bdy_getmoms(d,it,zf,rmcb,ymcb,nmom)
        use trdatbuf_obj_tsc
         implicit NONE
        type (trdatbuf) :: d
        integer, intent(in) :: it     ! time bin index
        real*8, intent(in) :: zf      ! time interpolation factor
        integer, intent(in) :: nmom      ! number of moments
        real*8, dimension(0:nmom,2) :: rmcb,ymcb  ! boundary harmonics 0:nmom, cos:sin
      end subroutine tdb_bdy_getmoms

  end interface
end module trdatbuf_iface
