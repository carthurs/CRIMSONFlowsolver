!$$$	subroutine e3DC (g1yi,   g2yi,   g3yi,  
!$$$     &                   giju,   tauM,    A0,     raLS,
!$$$     &			 rtLS,   giju,   DC,     ri,
!$$$     &                   rmi,    stiff, A0DC)
!$$$c
!$$$c----------------------------------------------------------------------
!$$$c
!$$$c This routine calculates the contribution of the Discontinuity-
!$$$c Capturing operator to RHS and preconditioner.
!$$$c
!$$$c  g1yi   (nflow,npro)           : grad-y in direction 1
!$$$c  g2yi   (nflow,npro)           : grad-y in direction 2
!$$$c  g3yi   (nflow,npro)           : grad-y in direction 3
!$$$c  A0     (nsymdf,npro)          : A0 matrix (Symm. storage)
!$$$c  raLS   (npro)                 : square of LS residual (A0inv norm)
!$$$c  rtLS   (npro)                 : square of LS residual (Tau norm)
!$$$c  giju    (6,npro)              : metric matrix
!$$$c  DC     (ngauss,npro)          : discontinuity-capturing factor
!$$$c  intp				 : integration point number
!$$$c
!$$$c output:
!$$$c  ri     (nflow*(nsd+1),npro)   : partial residual
!$$$c  rmi    (nflow*(nsd+1),npro)   : partial modified residual
!$$$c  stiff  (nsymdf,6,npro)       : diffusivity matrix
!$$$c  DC     (npro)                : discontinuity-capturing factor
!$$$c
!$$$c
!$$$c Zdenek Johan, Summer 1990. (Modified from e2dc.f)
!$$$c Zdenek Johan, Winter 1991. (Recoded)
!$$$c Zdenek Johan, Winter 1991. (Fortran 90)
!$$$c----------------------------------------------------------------------
!$$$c
!$$$	use phcommonvars
!$$$c
!$$$        dimension g1yi(npro,nflow),          g2yi(npro,nflow),
!$$$     &            g3yi(npro,nflow),          A0(npro,5,5),
!$$$     &            raLS(npro),                rtLS(npro),
!$$$     &            giju(npro,6),              DC(npro,ngauss),
!$$$     &            ri(npro,nflow*(nsd+1)),    rmi(npro,nflow*(nsd+1)),
!$$$     &            stiff(npro,3*nflow,3*nflow),itmp(npro)
!$$$c
!$$$
!$$$        dimension ggyi(npro,nflow),         gAgyi(npro,15),
!$$$     &            gnorm(npro),              A0gyi(npro,15),
!$$$     &            yiA0DCyj(npro,6),         A0DC(npro,4)
!$$$c
!$$$c ... -----------------------> initialize <----------------------------
!$$$c
!$$$        A0gyi    = zero
!$$$        gAgyi    = zero
!$$$        yiA0DCyj = zero
!$$$        DC       = zero
!$$$c.... ----------------------->  global gradient  <----------------------
!$$$c
!$$$c.... calculate (A0 y_,j) --> A0gyi
!$$$c
!$$$c  A0 Y_{,1}
!$$$c
!$$$        A0gyi( :,1) = A0(:,1,1)*g1yi(:,1)
!$$$     &              + A0(:,1,2)*g1yi(:,2)
!$$$     &              + A0(:,1,3)*g1yi(:,3)
!$$$     &              + A0(:,1,4)*g1yi(:,4)
!$$$     &              + A0(:,1,5)*g1yi(:,5)      
!$$$        A0gyi( :,2) = A0(:,2,1)*g1yi(:,1)
!$$$     &              + A0(:,2,2)*g1yi(:,2)
!$$$     &              + A0(:,2,3)*g1yi(:,3)
!$$$     &              + A0(:,2,4)*g1yi(:,4)
!$$$     &              + A0(:,2,5)*g1yi(:,5)
!$$$        A0gyi( :,3) = A0(:,3,1)*g1yi(:,1)
!$$$     &              + A0(:,3,2)*g1yi(:,2)
!$$$     &              + A0(:,3,3)*g1yi(:,3)
!$$$     &              + A0(:,3,4)*g1yi(:,4)
!$$$     &              + A0(:,3,5)*g1yi(:,5)
!$$$        A0gyi( :,4) = A0(:,4,1)*g1yi(:,1)
!$$$     &              + A0(:,4,2)*g1yi(:,2)
!$$$     &              + A0(:,4,3)*g1yi(:,3)
!$$$     &              + A0(:,4,4)*g1yi(:,4)
!$$$     &              + A0(:,4,5)*g1yi(:,5)
!$$$        A0gyi( :,5) = A0(:,5,1)*g1yi(:,1)
!$$$     &              + A0(:,5,2)*g1yi(:,2)
!$$$     &              + A0(:,5,3)*g1yi(:,3)
!$$$     &              + A0(:,5,4)*g1yi(:,4)
!$$$     &              + A0(:,5,5)*g1yi(:,5)
!$$$c
!$$$c  A0 Y_{,2}
!$$$c
!$$$        A0gyi( :,6) = A0(:,1,1)*g2yi(:,1)
!$$$     &              + A0(:,1,2)*g2yi(:,2)
!$$$     &              + A0(:,1,3)*g2yi(:,3)
!$$$     &              + A0(:,1,4)*g2yi(:,4)
!$$$     &              + A0(:,1,5)*g2yi(:,5)
!$$$        A0gyi( :,7) = A0(:,2,1)*g2yi(:,1)
!$$$     &              + A0(:,2,2)*g2yi(:,2)
!$$$     &              + A0(:,2,3)*g2yi(:,3)
!$$$     &              + A0(:,2,4)*g2yi(:,4)
!$$$     &              + A0(:,2,5)*g2yi(:,5)
!$$$        A0gyi( :,8) = A0(:,3,1)*g2yi(:,1)
!$$$     &              + A0(:,3,2)*g2yi(:,2)
!$$$     &              + A0(:,3,3)*g2yi(:,3)
!$$$     &              + A0(:,3,4)*g2yi(:,4)
!$$$     &              + A0(:,3,5)*g2yi(:,5)
!$$$        A0gyi( :,9) = A0(:,4,1)*g2yi(:,1)
!$$$     &              + A0(:,4,2)*g2yi(:,2)
!$$$     &              + A0(:,4,3)*g2yi(:,3)
!$$$     &              + A0(:,4,4)*g2yi(:,4)
!$$$     &              + A0(:,4,5)*g2yi(:,5)
!$$$        A0gyi(:,10) = A0(:,5,1)*g2yi(:,1)
!$$$     &              + A0(:,5,2)*g2yi(:,2)
!$$$     &              + A0(:,5,3)*g2yi(:,3)
!$$$     &              + A0(:,5,4)*g2yi(:,4)
!$$$     &              + A0(:,5,5)*g2yi(:,5)
!$$$c
!$$$c  A0 Y_{,3}
!$$$c
!$$$        A0gyi(:,11) = A0(:,1,1)*g3yi(:,1)
!$$$     &              + A0(:,1,2)*g3yi(:,2)
!$$$     &              + A0(:,1,3)*g3yi(:,3)
!$$$     &              + A0(:,1,4)*g3yi(:,4)
!$$$     &              + A0(:,1,5)*g3yi(:,5)
!$$$        A0gyi(:,12) = A0(:,2,1)*g3yi(:,1)
!$$$     &              + A0(:,2,2)*g3yi(:,2)
!$$$     &              + A0(:,2,3)*g3yi(:,3)
!$$$     &              + A0(:,2,4)*g3yi(:,4)
!$$$     &              + A0(:,2,5)*g3yi(:,5)
!$$$        A0gyi(:,13) = A0(:,3,1)*g3yi(:,1)
!$$$     &              + A0(:,3,2)*g3yi(:,2)
!$$$     &              + A0(:,3,3)*g3yi(:,3)
!$$$     &              + A0(:,3,4)*g3yi(:,4)
!$$$     &              + A0(:,3,5)*g3yi(:,5)
!$$$        A0gyi(:,14) = A0(:,4,1)*g3yi(:,1)
!$$$     &              + A0(:,4,2)*g3yi(:,2)
!$$$     &              + A0(:,4,3)*g3yi(:,3)
!$$$     &              + A0(:,4,4)*g3yi(:,4)
!$$$     &              + A0(:,4,5)*g3yi(:,5)
!$$$        A0gyi(:,15) = A0(:,5,1)*g3yi(:,1)
!$$$     &              + A0(:,5,2)*g3yi(:,2)
!$$$     &              + A0(:,5,3)*g3yi(:,3)
!$$$     &              + A0(:,5,4)*g3yi(:,4)
!$$$     &              + A0(:,5,5)*g3yi(:,5)
!$$$c
!$$$c.... calculate (giju A0 y_,j) --> gAgyi
!$$$c
!$$$
!$$$        gAgyi( :,1) = giju(:,1)*A0gyi( :,1)
!$$$     &              + giju(:,4)*A0gyi( :,6)
!$$$     &              + giju(:,5)*A0gyi(:,11)
!$$$
!$$$        gAgyi( :,2) = giju(:,1)*A0gyi( :,2)
!$$$     &              + giju(:,4)*A0gyi( :,7)
!$$$     &              + giju(:,5)*A0gyi(:,12)
!$$$
!$$$	gAgyi( :,3) = giju(:,1)*A0gyi( :,3)
!$$$     &              + giju(:,4)*A0gyi( :,8)
!$$$     &              + giju(:,5)*A0gyi(:,13)
!$$$
!$$$	gAgyi( :,4) = giju(:,1)*A0gyi( :,4)
!$$$     &              + giju(:,4)*A0gyi( :,9)
!$$$     &              + giju(:,5)*A0gyi(:,14)
!$$$
!$$$	gAgyi( :,5) = giju(:,1)*A0gyi( :,5)
!$$$     &              + giju(:,4)*A0gyi(:,10)
!$$$     &              + giju(:,5)*A0gyi(:,15)
!$$$
!$$$	gAgyi( :,6) = giju(:,4)*A0gyi( :,1)
!$$$     &              + giju(:,2)*A0gyi( :,6)
!$$$     &              + giju(:,6)*A0gyi(:,11)
!$$$
!$$$	gAgyi( :,7) = giju(:,4)*A0gyi( :,2)
!$$$     &              + giju(:,2)*A0gyi( :,7)
!$$$     &              + giju(:,6)*A0gyi(:,12)
!$$$
!$$$	gAgyi( :,8) = giju(:,4)*A0gyi( :,3)
!$$$     &              + giju(:,2)*A0gyi( :,8)
!$$$     &              + giju(:,6)*A0gyi(:,13)
!$$$
!$$$	gAgyi( :,9) = giju(:,4)*A0gyi( :,4)
!$$$     &              + giju(:,2)*A0gyi( :,9)
!$$$     &              + giju(:,6)*A0gyi(:,14)
!$$$
!$$$	gAgyi(:,10) = giju(:,4)*A0gyi( :,5)
!$$$     &              + giju(:,2)*A0gyi(:,10)
!$$$     &              + giju(:,6)*A0gyi(:,15)
!$$$
!$$$	gAgyi(:,11) = giju(:,5)*A0gyi( :,1)
!$$$     &              + giju(:,6)*A0gyi( :,6)
!$$$     &              + giju(:,3)*A0gyi(:,11)
!$$$
!$$$	gAgyi(:,12) = giju(:,5)*A0gyi( :,2)
!$$$     &              + giju(:,6)*A0gyi( :,7)
!$$$     &              + giju(:,3)*A0gyi(:,12)
!$$$
!$$$	gAgyi(:,13) = giju(:,5)*A0gyi( :,3)
!$$$     &              + giju(:,6)*A0gyi( :,8)
!$$$     &              + giju(:,3)*A0gyi(:,13)
!$$$
!$$$	gAgyi(:,14) = giju(:,5)*A0gyi( :,4)
!$$$     &              + giju(:,6)*A0gyi( :,9)
!$$$     &              + giju(:,3)*A0gyi(:,14)
!$$$
!$$$	gAgyi(:,15) = giju(:,5)*A0gyi( :,5)
!$$$     &              + giju(:,6)*A0gyi(:,10)
!$$$     &              + giju(:,3)*A0gyi(:,15)
!$$$c	
!$$$c... the denominator term of the DC factor
!$$$c... evaluation of the term  Y,i.A0DC Y,j 
!$$$c
!$$$        yiA0DCyj(:,1) = A0DC(:,1)*g1yi(:,1)**2
!$$$     &                + two*g1yi(:,1)*A0DC(:,2)*g1yi(:,5)
!$$$     &                + A0DC(:,3)*g1yi(:,2)**2
!$$$     &                + A0DC(:,3)*g1yi(:,3)**2
!$$$     &                + A0DC(:,3)*g1yi(:,4)**2
!$$$     &                + A0DC(:,4)*g1yi(:,5)**2
!$$$
!$$$        yiA0DCyj(:,2) = A0DC(:,1)*g2yi(:,1)**2
!$$$     &                + two*g2yi(:,1)*A0DC(:,2)*g2yi(:,5)
!$$$     &                + A0DC(:,3)*g2yi(:,2)**2
!$$$     &                + A0DC(:,3)*g2yi(:,3)**2
!$$$     &                + A0DC(:,3)*g2yi(:,4)**2
!$$$     &                + A0DC(:,4)*g2yi(:,5)**2
!$$$
!$$$        yiA0DCyj(:,3) = A0DC(:,1)*g3yi(:,1)**2
!$$$     &                + two*g3yi(:,1)*A0DC(:,2)*g3yi(:,5)
!$$$     &                + A0DC(:,3)*g3yi(:,2)**2
!$$$     &                + A0DC(:,3)*g3yi(:,3)**2
!$$$     &                + A0DC(:,3)*g3yi(:,4)**2
!$$$     &                + A0DC(:,4)*g3yi(:,5)**2
!$$$
!$$$        yiA0DCyj(:,4) = g1yi(:,1)*A0DC(:,1)*g2yi(:,1)
!$$$     &                + g1yi(:,1)*A0DC(:,2)*g2yi(:,5)
!$$$     &                + g1yi(:,2)*A0DC(:,3)*g2yi(:,2)
!$$$     &                + g1yi(:,3)*A0DC(:,3)*g2yi(:,3)
!$$$     &                + g1yi(:,4)*A0DC(:,3)*g2yi(:,4)
!$$$     &                + g1yi(:,5)*A0DC(:,2)*g2yi(:,1)
!$$$     &                + g1yi(:,5)*A0DC(:,4)*g2yi(:,5)
!$$$
!$$$        yiA0DCyj(:,5) = g1yi(:,1)*A0DC(:,1)*g3yi(:,1)
!$$$     &                + g1yi(:,1)*A0DC(:,2)*g3yi(:,5)
!$$$     &                + g1yi(:,2)*A0DC(:,3)*g3yi(:,2)
!$$$     &                + g1yi(:,3)*A0DC(:,3)*g3yi(:,3)
!$$$     &                + g1yi(:,4)*A0DC(:,3)*g3yi(:,4)
!$$$     &                + g1yi(:,5)*A0DC(:,2)*g3yi(:,1)
!$$$     &                + g1yi(:,5)*A0DC(:,4)*g3yi(:,5)
!$$$
!$$$        yiA0DCyj(:,6) = g2yi(:,1)*A0DC(:,1)*g3yi(:,1)
!$$$     &                + g2yi(:,1)*A0DC(:,2)*g3yi(:,5)
!$$$     &                + g2yi(:,2)*A0DC(:,3)*g3yi(:,2)
!$$$     &                + g2yi(:,3)*A0DC(:,3)*g3yi(:,3)
!$$$     &                + g2yi(:,4)*A0DC(:,3)*g3yi(:,4)
!$$$     &                + g2yi(:,5)*A0DC(:,2)*g3yi(:,1)
!$$$     &                + g2yi(:,5)*A0DC(:,4)*g3yi(:,5)
!$$$c
!$$$c.... ------------------------->  DC factor  <--------------------------
!$$$c
!$$$	if ((ires .ne. 2) .or. (Jactyp .eq. 1)) then
!$$$c
!$$$c.... calculate 2-norm of Grad-local-V with respect to A0
!$$$c
!$$$c.... DC-mallet
!$$$c
!$$$	  if (iDC .eq. 1) then
!$$$c
!$$$	    fact = one
!$$$	    if (ipord .eq. 2)  fact = 0.9
!$$$	    if (ipord .eq. 3) fact = 0.75
!$$$	
!$$$c
!$$$            gnorm = one / (
!$$$     &              giju(:,1)*yiA0DCyj(:,1)
!$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
!$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
!$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
!$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
!$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
!$$$     &            + epsM  )
!$$$c
!$$$	    DC(:,intp)=dim((fact*sqrt(raLS*gnorm)),(rtLS*gnorm))
!$$$c
!$$$c.... flop count
!$$$c
!$$$	    flops = flops + 46*npro
!$$$c
!$$$	  endif
!$$$c
!$$$c.... DC-quadratic
!$$$c
!$$$	  if (iDC .eq. 2) then
!$$$c
!$$$            gnorm = one / (
!$$$     &              giju(:,1)*yiA0DCyj(:,1)
!$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
!$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
!$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
!$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
!$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
!$$$     &            + epsM  )
!$$$         
!$$$c
!$$$	    DC(:,intp) = two * rtLS * gnorm
!$$$c
!$$$c.... flop count
!$$$c
!$$$	    flops = flops + 36*npro
!$$$c
!$$$	  endif
!$$$c
!$$$c.... DC-min
!$$$c
!$$$	  if (iDC .eq. 3) then
!$$$c
!$$$	    fact = one
!$$$	    if (ipord .eq. 2)  fact = pt5
!$$$c
!$$$            gnorm = one / (
!$$$     &              giju(:,1)*yiA0DCyj(:,1)
!$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
!$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
!$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
!$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
!$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
!$$$     &            + epsM  )
!$$$
!$$$c
!$$$	    DC(:,intp) = min( dim(fact * sqrt(raLS * gnorm),
!$$$     &                       rtLS * gnorm), two * rtLS * gnorm )
!$$$c
!$$$c.... flop count
!$$$c
!$$$	    flops = flops + 48*npro
!$$$c
!$$$	  endif
!$$$c
!$$$	endif
!$$$c
!$$$c.... ---------------------------->  RHS  <----------------------------
!$$$c
!$$$c.... add the contribution of DC to ri and/or rmi
!$$$c
!$$$c.... ires = 1 or 3
!$$$c
!$$$	if ((ires .eq. 1) .or. (ires .eq. 3)) then
!$$$c
!$$$	  ri ( :,1) = ri ( :,1) + DC(:,intp) * gAgyi( :,1)
!$$$	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
!$$$	  ri ( :,2) = ri ( :,2) + DC(:,intp) * gAgyi( :,2)
!$$$	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
!$$$	  ri ( :,3) = ri ( :,3) + DC(:,intp) * gAgyi( :,3)
!$$$	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
!$$$	  ri ( :,4) = ri ( :,4) + DC(:,intp) * gAgyi( :,4)
!$$$	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
!$$$	  ri ( :,5) = ri ( :,5) + DC(:,intp) * gAgyi( :,5)
!$$$	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
!$$$c
!$$$	  ri ( :,6) = ri ( :,6) + DC(:,intp) * gAgyi( :,6)
!$$$	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
!$$$	  ri ( :,7) = ri ( :,7) + DC(:,intp) * gAgyi( :,7)
!$$$	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
!$$$	  ri ( :,8) = ri ( :,8) + DC(:,intp) * gAgyi( :,8)
!$$$	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
!$$$	  ri ( :,9) = ri ( :,9) + DC(:,intp) * gAgyi( :,9)
!$$$	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
!$$$	  ri (:,10) = ri (:,10) + DC(:,intp) * gAgyi(:,10)
!$$$	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
!$$$c
!$$$	  ri (:,11) = ri (:,11) + DC(:,intp) * gAgyi(:,11)
!$$$	  rmi(:,11) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
!$$$	  ri (:,12) = ri (:,12) + DC(:,intp) * gAgyi(:,12)
!$$$	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
!$$$	  ri (:,13) = ri (:,13) + DC(:,intp) * gAgyi(:,13)
!$$$	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
!$$$	  ri (:,14) = ri (:,14) + DC(:,intp) * gAgyi(:,14)
!$$$	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
!$$$	  ri (:,15) = ri (:,15) + DC(:,intp) * gAgyi(:,15)
!$$$	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
!$$$c
!$$$	  flops = flops + 45*npro
!$$$c
!$$$	endif
!$$$c
!$$$c.... ires = 2
!$$$c
!$$$	if (ires .eq. 2) then
!$$$c
!$$$	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
!$$$	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
!$$$	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
!$$$	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
!$$$	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
!$$$c
!$$$	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
!$$$	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
!$$$	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
!$$$	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
!$$$	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
!$$$c
!$$$	  rmi(:,11) = rmi(:,11) + DC(:,intp) * gAgyi(:,11)
!$$$	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
!$$$	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
!$$$	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
!$$$	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
!$$$c
!$$$	  flops = flops + 30*npro
!$$$c
!$$$	endif
!$$$c
!$$$c.... ------------------------->  Stiffness  <--------------------------
!$$$c
!$$$c.... add the contribution of DC to stiff
!$$$c
!$$$	if (iprec .eq. 1) then
!$$$	     nflow2=two*nflow
!$$$       do j = 1, nflow
!$$$          do i = 1, nflow
!$$$             itmp(:)=A0(:,i,j)*DC(:,intp)
!$$$c
!$$$c.... add (DC g^1 A0) to stiff [1,1]
!$$$c
!$$$             stiff(:,i,j) = stiff(:,i,j) 
!$$$     &                    + itmp(:)*giju(:,1)
!$$$c
!$$$c.... add (DC g^1 A0) to stiff [1,2]
!$$$c
!$$$
!$$$             stiff(:,i,j+nflow) = stiff(:,i,j+nflow) 
!$$$     &                    + itmp(:)*giju(:,4)
!$$$c
!$$$c.... add (DC g^1 A0) to stiff [1,3]
!$$$c
!$$$
!$$$             stiff(:,i,j+nflow2) = stiff(:,i,j+nflow2) 
!$$$     &                    + itmp(:)*giju(:,5)
!$$$
!$$$c.... add (DC g^1 A0) to stiff [2,1] (similarly below)
!$$$c
!$$$
!$$$             stiff(:,i+nflow,j) = stiff(:,i+nflow,j) 
!$$$     &                    + itmp(:)*giju(:,4)
!$$$
!$$$             stiff(:,i+nflow,j+nflow) = stiff(:,i+nflow,j+nflow) 
!$$$     &                    + itmp(:)*giju(:,2)
!$$$
!$$$             stiff(:,i+nflow,j+nflow2) = stiff(:,i+nflow,j+nflow2) 
!$$$     &                    + itmp(:)*giju(:,6)
!$$$
!$$$             stiff(:,i+nflow2,j) = stiff(:,i+nflow2,j) 
!$$$     &                    + itmp(:)*giju(:,5)
!$$$
!$$$             stiff(:,i+nflow2,j+nflow) = stiff(:,i+nflow2,j+nflow) 
!$$$     &                    + itmp(:)*giju(:,6)
!$$$
!$$$             stiff(:,i+nflow2,j+nflow2) = stiff(:,i+nflow2,j+nflow2) 
!$$$     &                    + itmp(:)*giju(:,3)
!$$$          enddo
!$$$       enddo
!$$$c
!$$$c.... flop count
!$$$c
!$$$	  flops = flops + 210*npro
!$$$c
!$$$c.... end of stiffness
!$$$c
!$$$	endif
!$$$c
!$$$c.... return
!$$$c
!$$$	return
!$$$	end
!$$$c
!
        subroutine e3dcSclr ( gradS,    giju,     gGradS, & 
                              rLS,      tauS,     srcR, & 
                              dcFct)
!
!
!----------------------------------------------------------------------
!
! This routine calculates the contribution of the Discontinuity-
! Capturing operator to RHS and preconditioner for the scalar solve.
!
!  g1yti   (nflow,npro)           : grad-y in direction 1
!  g2yti   (nflow,npro)           : grad-y in direction 2
!  g3yti   (nflow,npro)           : grad-y in direction 3
!  A0     (nsymdf,npro)          : A0 matrix (Symm. storage)
!  raLS   (npro)                 : square of LS residual (A0inv norm)
!  rtLS   (npro)                 : square of LS residual (Tau norm)
!  giju    (6,npro)              : metric matrix
!  DC     (ngauss,npro)          : discontinuity-capturing factor
!  intp				 : integration point number
!
! output:
!  ri     (nflow*(nsd+1),npro)   : partial residual
!  rmi    (nflow*(nsd+1),npro)   : partial modified residual
!  stiff  (nsymdf,6,npro)       : diffusivity matrix
!  DC     (npro)                : discontinuity-capturing factor
!
!
! Zdenek Johan, Summer 1990. (Modified from e2dc.f)
! Zdenek Johan, Winter 1991. (Recoded)
! Zdenek Johan, Winter 1991. (Fortran 90)
!----------------------------------------------------------------------
!
	 use phcommonvars  
 	IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension gradS(npro,nsd),            gGradS(npro,nsd), & 
                  rLS(npro),                  tauS(npro), & 
                  giju(npro,6),               dcFct(npro), & 
                  srcR(npro)
!
!.... Form GijUp gradS and  gradS . GijUp gradS (store in dcFct)
!
	
	    gGradS(:,1) = GijU(:,1) * gradS(:,1) &
     			+ GijU(:,4) * gradS(:,2) &
     			+ GijU(:,6) * gradS(:,3)
	    gGradS(:,2) = GijU(:,4) * gradS(:,1) &
     			+ GijU(:,2) * gradS(:,2) &
     			+ GijU(:,5) * gradS(:,3)
	    gGradS(:,3) = GijU(:,6) * gradS(:,1) &
     			+ GijU(:,5) * gradS(:,2) &
     			+ GijU(:,3) * gradS(:,3)
!
	    dcFct(:)    = gradS(:,1) * gGradS(:,1) &
     		        + gradS(:,2) * gGradS(:,2) &
     		        + gradS(:,3) * gGradS(:,3) &
     		        + epsM
	
	    dcFct(:) = 1.0/ dcFct(:)
!
!.... Form pdeRes 2-norm / gradT 2-norm
!

	    dcFct  = dcFct * (rLS - srcR) ** 2 
!
!.... ------------------------->  DC factor  <------------------------
!
!.... DC-mallet
!
	    if (idcsclr(1) .eq. 1) then
!       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
	       if (ipord .eq. 3) fact = 0.75
!       
!$$$  dcFct(:)=dim((fact*sqrt(dcFct(:))),(tauS(:)*dcFct(:))) !not work
                                                          !with all compilers
	       dcFct(:)=max(0.0,(fact*sqrt(dcFct(:)))-(tauS(:)*dcFct(:)))
!
	    endif
!       
!       
!....   DC-quadratic
!       
	    if (idcsclr(1) .eq. 2) then
!       
	       dcFct(:) = two * tauS(:) * dcFct(:)
!       
	    endif
!       
!....   DC-min
!       
	    if (idcsclr(1) .eq. 3) then
!       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
!       
          dcFct(:) = min( max(0.0, (fact * sqrt(dcFct(:)) - & 
      	             tauS(:)*dcFct(:)) ), two * tauS(:) * dcFct(:))
!       
	    endif
!
!.... Scale the gGradT for residual formation
!	
	    gGradS(:,1) = dcFct(:) * gGradS(:,1)
	    gGradS(:,2) = dcFct(:) * gGradS(:,2)
	    gGradS(:,3) = dcFct(:) * gGradS(:,3)
	


	return
	end
!
