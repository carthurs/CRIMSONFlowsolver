      subroutine scatnu (ien, strl, xmudmi, xnut, shp)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension  ien(npro,nshl),       strl(npro,ngauss), &
                 xmudmi(npro,ngauss),       shp(nshl,ngauss)
      dimension  xnut(numnp)

      xmudmi=zero

      if(iLES.eq.5) return  ! Debugging with zero-ed model

      do in = 1,nshl
      do int = 1, ngauss
        xmudmi(:,int) = xmudmi(:,int) + xnut(ien(:,in)) * strl(:,int) &
              *shp(in,int)
      enddo  
      enddo
!
!  local clipping
!
      rmu=datmat(1,2,1)
      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
      xmudmi=max(xmudmi, zero) ! don't let (xmudmi) < 0
!      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
      return
      end
