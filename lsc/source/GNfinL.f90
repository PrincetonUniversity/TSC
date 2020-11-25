!
!----------------------------------------------------------------------
!
      SUBROUTINE GNfinL
      USE gnuI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      if ( level .ge. 2 ) then
           write (iUnit,'(a)') '# empty graph trick, start #'
           write (iUnit, 600)
 600       format('set origin 0.80,0.0',/,'set size    0.15, 0.15',/,    &  
     &       'set xrange [-1:1]',/,      'set yrange [ 0.:1]',/,         &  
     &       'set nozeroaxis',/,         'set noborder',/,               &  
     &       'set noxtics',/,            'set noytics',/,                &  
     &       'plot ''-'' with impulses',/,                               &  
     &       '  0 0',/,                                                  &  
     &       'e',/,                                                      &  
     &       'set nomultiplot',/,                                        &  
     &       'reset')
         write (iUnit,'(a)') 'set nokey'
         write (iUnit,'(a)') 'set title'
           write (iUnit,'(a)') '# empty graph trick, stop #'
           level=1
      else
           write (iUnit,'(a)') '# end of file #'
           close (iUnit)
           level=0
      endif
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
