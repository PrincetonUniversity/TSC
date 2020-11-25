      include 'param.i'
      character*5 ident
      character*8 name
      common/only/ a(ndim),b(ndim),c(ndim),d(ndim),e(ndim),
     1 f(ndim),psi(ndim),r(ndim),aj0(ndim),ajcheck(ndim),
     2  eta(ndim),ajcd(ndim),etai(ndim),bpsq(ndim),
     3  anre(ndim),ane(ndim),etal(ndim),ajb(ndim),aout(ndim)
      common/scalar/ dr
      common/iscalar/ iwkid,iwtype,lunit,ierrf,ident(3),maxvar,
     1   name(20)
c
