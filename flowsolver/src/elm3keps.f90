!C===========================================================================
!C
!C "elm3Keps.f":  3D k-epsilon (LB model) turbulence element formation.
!C
!C Original: Farzin Shakib (Nov 98)
!C K-epsilon modification by Je-Hoon, Kim (Jan 98)
!C Wlicox's book, pp.140, Lam-Bremhorst model
!C
!C Includes k-epsilons stagnation abnormally correction (Feb 99)
!C ( Durbin, P. A., "On the k-3 stagnation point
!C   abnormally", pp.89-90, 1996, Int. J. of Heat and Fluid flow).
!C


!... Note

!..  1. K-epsilon model seems more sensitive to grid resolution than one equation models.
!        It seems good (fine enough) grid seldom get blow-up ('nan' error) while coarse grids
!        sometimes blow up. It is strongly recommended to have very fine grid near the wall
!        at least y+ (1st grid off the wall) less than or = 0.5. 

!..  2. when it blew up while converging to the solution, alternative srcJac's were
!..     used to prevent blow up. It eliminated blow up problem, but it does not seem
!..     to be a good Jacobian to solve the problem from the very beginning (just good
!..     as a restarter).
!..     Further numerical test will be done to find a better Jacobian
!..     which is good for either case, I have not finished it yet.
!..     I have coded all possible terms to play with and am studying the behavior now. 

!..  3. Durbin stated that we could limit the T2nd (second time scale) to control
!..     the abnormal growth of k , but it does not seem to help much. 

!..  4. About good IC's for K-epsilon model,
!..     K_(O)      =   k+ max(=4.5~6) * (u_tau)^2
!..     eps_(O)    = eps+ max(=0.1~0.25) * (u_tau)^4 / kinematic_viscosity
!..     velocity(0)= u_tau *( 1/0.418*Log(y+ max)+5.5 ) or mass velocity if it can be estimated.
!..     mu_t(0)    = (500~1000) * mu



!..     (example) for periodic channel problem 

!..             when  u_tau=0.707,  k+   max = 6,    k IC   =    3
!..                                 eps+ max = 0.1,  eps IC = 1250 
!..                                                  vel IC ~   20

!..             when  u_tau=1,       k+   max = 6,    k IC  =    6
!..                                 eps+ max = 0.1,  eps IC = 5000 
!..                                                  vel IC ~   30
!..      u_tau=1's IC's works for lower dpdx cases.


!..     I used possible maximum for k and epsilon.
!..     From one of review papers on k-epsilon models, we know approximate k+ max
!..     and epslon+ max for attached flows. 
!..     Shear velocity (utau) could also be estimated for particular probelm which
!..     actually set k and epsilon along with another flow properties.
!..     Too small epsilon or too small mu_t might let k grow unchecked therefore
!..     make the flow field unstable.

!..  5. Restart works better from higher Reynolds number to lower Reynolds number flows but
!..     not the other way (Before Durbin's modification, and I haven't check after Durbin's
!..     modification).

!..  6. In terms of terminology and naming convention, I pretty much followed DR. Shakib's
!..     way. If there are any confusion, pls let me know.


!..  7. Slight modifications are done to Jacobians (basically sign control) to accelerate 
!..     convergence and to prevent a blow-up before the solutions converge


!..  8. Periodic channel problem tells :
!..     The lower bound of y+ of the 1st node is 0.6 and 12 nodes within y+ = 20 (near wall region)
!..     (STAR-CD recommend y+~1 and 20 nodes within y+=20).
!..     Our 1st y+ requirement is sligtly more than what STAR-CD recommends. 


!..  9. Fin tube problem (pls neglect the sample problem. it will not converge) tells :
!..     By inspecting y+ resolution, when the solver blowed up, there were
!..     about 6~7 points which had y+ (1st node) more than 0.6 near the stagnation points.
!..     The GR's also did not go down below GR of 1.e-4.

!..     After fixing y+ resolution near the stagnation points, it converged well (below 1.e-4).
!..
!..     (example) IC selection for fin tube
 
!..              To find out candidiate Initial Condition for k and epsilon, I assumed that
!..              the wall shear stress could reach upto 10 times the average wall shear stress
!..              at the stagnation point (I think it will vary of course, depending on the problem).
!..              
!..               u_tau_average can be estimated from experimaental data.
!..               u_tau_max = SQRT(10) * u_tau_average
!..               k+ max ~ 6, epsilon max ~0.25
!..               
!C===========================================================================


!============================================================================
!
! "fElm3KepsCoef":
!
!============================================================================
      subroutine elm3keps    (kay,	epsilon, &
           dwi,	gradV, &
           srcRat1,src1,	srcJac)
!
!.... Data declaration
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      real*8	kay(npro),		epsilon(npro), &
                dwi(npro)
      real*8	gradV(npro,3,3)

      real*8	rmu(npro), &
                rho(npro)
      real*8	srcRat(npro,2), &
      		src(npro,2),		srcJac(npro,4), &
                srcRat1(npro),          src1(npro)
!
      integer advdiff
      integer	e
      real*8  fct1
      real*8  tmp1,   tmp2  ,  tmp0
!
      real*8  k,    kInv,   kq,     eps,   y,    epsInv,   epsP2, &
                ss,   mut, mut_k, mut_eps,  rat,      jac
    

      real*8  Ce1,   Ce2,  epsP2Inv,   CmuKE,  sigmaKE,       nu,  &
                kk,  kqInv,     nuInv,     Rey,     Ret,    Rey_k,   &
                Ret_k,        Ret_eps,     ff1,      f2,    kkInv, &
                RetInv,       fmukeInv,  fmukeP2Inv

      real*8  fmuke,  fmuke_k,  fmuke_eps

!...    Addings due to Durbin's correction 

      real*8  T1st, T2nd, T2ndInv, two3rdCmuInv, tri8thq, &
                ssq, ssqInv, two3rdq, pk1, pk2,  &
                pk, pk_k, pk_eps, Tscale,Tinv, &
                Tscale_k, Tscale_eps, ff1_k, ff1_eps, f2_k,f2_eps, &
                ff1Inv, f2Inv

 

!
!	include "elmf.h"
!
!.... Compute src and its jacobians
!
!.... set Lam-Bremhorst' k-epsilon Model constants 
      rho(:)=datmat(1,1,1)
      rmu(:)=datmat(1,2,1)
      Ce1             = 1.44
      Ce2             = 1.92
      CmuKE           = 0.09
      sigmaKE         = 1.3
      
      two3rdCmuInv    = 2./3./CmuKE
      tri8thq         = SQRT(3./8.)
      two3rdq         = SQRT(2./3.)
!
      advdiff = 0
      if(advdiff.eq.0)then ! not advection-diffusion
      do e = 1, npro
         
         nuInv           = rho(e)/rmu(e)
         k		    = abs(kay(e))
         if (k.lt.1.e-32) k=0
         eps		    = abs(epsilon(e))
         y		    = dwi(e)

!--------------patch
!            if(k.gt.0.45) k=0.45
!            if(eps.gt.2158) eps=2158
!--------------------------------

!
         kInv	    = 0
         kq	            = sqrt(k)
         kqInv           = 0
         kkInv           = 0
 
!.... limiting k.  Instead of saying k.ne.0

         if ( k .gt.1.e-32 ) then 
            kInv         = 1. / k
            kqInv        = 1./sqrt(k)
            kkInv        = kInv*kInv
         endif
            
         kk              =   k * k

         epsP2           = eps * eps
         epsInv	    = 0
         epsP2Inv        = 0


!.....  limiting epsilon.  Instead of saying eps.ne.0
             
         if ( eps .gt.1.e-32) then
            epsInv        = 1. / eps
            epsP2Inv      = epsInv*epsInv
         endif

!
         ss		= gradV(e,1,1) ** 2 &
      			+ gradV(e,2,2) ** 2 &
      			+ gradV(e,3,3) ** 2 &
      			+ 0.5 &
      			* ( (gradV(e,2,3) + gradV(e,3,2)) ** 2 &
      			  + (gradV(e,3,1) + gradV(e,1,3)) ** 2 &
      			  + (gradV(e,1,2) + gradV(e,2,1)) ** 2 )
!
         ssq                 = sqrt(ss)   
         ssqInv              = 0

         if(ss.ne.0) ssqInv  = 1./sqrt(ss)
       
         Rey                 = kq *      y * nuInv       
         Ret                 = kk * epsInv * nuInv     
         RetInv              = 0

!...     limitng Ret so that it does not get 'nan' error

         if(Ret.lt.1.d100.AND.Ret.gt.zero) RetInv=1./Ret

         Rey_k    =  0.5 * y * kqInv  * nuInv   
         Ret_k    =  2.  * k * epsInv * nuInv    
         Ret_eps  = -kk  * epsP2Inv   * nuInv      

         tmp1     = exp(-0.0165*Rey)          

         fmuKE    = (1. -tmp1) ** 2 * (1.+20.5*RetInv) ! fmu of Lam-Bremhorst

         fmuKEInv = 0.0
         fmuKEP2Inv = 0.0
!....   limiting fmuKE.  fmuke max ~ 1, and it could get very small near the wall.
 
         if(fmuKe.gt.1e-32) then
            fmuKEInv = 1./fmuke
            fmuKEP2Inv = fmuKEInv*fmuKEInv
         endif

         fmuKE_k  = 0.033*(1.-tmp1)*(1.+20.5*RetInv)*Rey_k*tmp1 &
                      -20.5 *(1.-tmp1)**2 * Ret_k*RetInv**2
   
         fmuKE_eps= -20.5*(1.-tmp1)**2* Ret_eps*RetInv**2     
         
         ff1      =  1. + ( 0.05*fmuKEInv) ** 3 ! f1 as in Lam-Bremhorst
         f2       =  1. - exp(- Ret ** 2 ) ! f2 as in Lam-Bremhorst     

         ff1Inv=zero
         f2Inv=zero
         if(ff1.gt.1.0e-32)    ff1Inv=1./ff1
         if(f2.gt.1.0e-32)     f2Inv =1./f2


         ff1_k    = -0.000375*fmuKE_k  *fmuKEInv**4
         ff1_eps  = -0.000375*fmuKE_eps*fmuKEInv**4
         
         f2_k     = 2.* Ret * Ret_k   * exp(-Ret**2)
         f2_eps   = 2.* Ret * Ret_eps * exp(-Ret**2)
         
         T1st = k * epsInv      ! 1st time scale
         T2nd = fmuKEInv*two3rdCmuInv*tri8thq*ssqInv ! 2nd time scale
         
!...     Depending on the choice of T (to limit k growth near stagnation)    

         if (T1st.lt.T2nd) then

            Tscale      = T1st
            Tinv        =  0

!                 if(Tscale.ne.0)  Tinv=1./Tscale

!...       Limiting time scale s.t. Tinv**4 not go over 1.e160

            if(Tscale.gt.1.e-32)  Tinv=1./Tscale
            Tscale_k= epsInv
            Tscale_eps= -k*epsP2Inv

         else

            Tscale      = T2nd
            Tinv        =  0

            if(Tscale.gt.1.e-32) Tinv=1./Tscale
            Tscale_k= two3rdCmuInv*tri8thq*ssqInv* &
                         (-fmuKEP2Inv*fmuKE_k) 
            Tscale_eps= two3rdCmuInv*tri8thq*ssqInv* &
                         (-fmuKEP2Inv*fmuKE_eps)  

         endif

!...      After Limiting all values, i feel it's unnecessary to limit jacobians
!...      since they are already limited.
!...      Two other routines which defines these quantities are the same




!...   recall that tmp1= exp(-0.0165*Rey) 


             
         mut  = rho(e)*CmuKE*fmuKE* k *Tscale

         mut_k = rho(e)*CmuKE*( &
                       fmuKE_k*k*Tscale &
                     + fmuKE    *Tscale &
                     + fmuKE  *k*Tscale_k &
                                  )

         mut_eps  = rho(e) * CmuKE *k*( &
              fmuKE_eps*Tscale+ &
              fmuKE    *Tscale_eps &
                                  ) 
        
         tmp0        = 2 * ss 
 

         pk1     = mut*tmp0
         pk2     = rho(e)*two3rdq*k*sqrt(ss)


         if (pk1.lt.pk2) then
            pk     = pk1
            pk_k   = mut_k  * tmp0 
            pk_eps = mut_eps * tmp0
         else
            pk   = pk2
            pk_k = kInv*pk2
            pk_eps =0
         endif


         src(e,1)    = pk    - rho(e) * eps      
         src(e,2)    = ff1 * Ce1 * pk * Tinv &
                         -rho(e)* Ce2 * f2 * eps * Tinv

         srcRat(e,1) = -kInv*(  pk    - rho(e) * eps   )   
         srcRat(e,2) = -epsInv*(   &
                          ff1 * Ce1 * pk * Tinv &
                         -rho(e)* Ce2 * f2 * eps * Tinv &
                               )

         srcJac(e,1) = -(pk_k)  ! jacobian for k PDE alone  
         srcJac(e,3) = -(pk_eps-rho(e)) ! d(Fsrck)/d(epsilon)
         srcJac(e,2) = -( &
                  ff1_k * Ce1 * pk    * Tinv &
                 +ff1   * Ce1 * pk_k  * Tinv &
                 -ff1   * Ce1 * pk    * Tscale_k * Tinv**2  &
                 -rho(e)  * Ce2 * f2_k  * eps * Tinv &
                 +rho(e)  * Ce2 * f2    * eps * Tscale_k * Tinv**2 ) ! do not touch  &
                            ! d(Fsrce)/d(k)
         srcJac(e,4) = -( &
                 ff1_eps * Ce1 * pk   * Tinv         &
                 +ff1     * Ce1 * pk_eps * Tinv &
              -ff1     * Ce1 * pk   * Tscale_eps * Tinv**2      &
                 -rho(e)    * Ce2 * f2_eps     * eps * Tinv   &
                 -rho(e)    * Ce2 * f2       * Tinv &
              +rho(e)    * Ce2 * f2* eps * Tscale_eps * Tinv**2   &
                    ) ! jacobian for epsilon PDE alone


       

      enddo
!
!.... Ensure positivity of srcJac
!
      do e = 1, npro


         srcJac(e,1)		= max( srcJac(e,1), srcRat(e,1), 0.d0 )
         srcJac(e,4)		= max( srcJac(e,4), srcRat(e,2), 0.d0 )

!            write(777,*) e,   srcJac(e,1),   srcJac(e,4)

!
         tmp1		= min( srcJac(e,1) * srcJac(e,4), &
      				       (srcJac(e,1)-srcRat(e,1)) * &
      				       (srcJac(e,4)-srcRat(e,2)) )
         tmp2		= srcJac(e,2) * srcJac(e,3)


         if ( tmp2 .gt. tmp1 ) then
            tmp2		= sqrt(tmp1/tmp2)
            srcJac(e,2)	= tmp2 * srcJac(e,2)
            srcJac(e,3)	= tmp2 * srcJac(e,3)
         endif
!	    srcJac(e,2)	= 0
!	    srcJac(e,3)	= 0
      enddo
      if(isclr.eq.1)then        ! kay
         srcrat1 = srcrat(:,1)
         src1 = src(:,1)
      else if (isclr.eq.2) then ! epsilon
         srcrat1 = srcrat(:,2)
         src1 = src(:,2)
      endif

      else ! Advection-diffusion
! Advection-diffusion case
      srcrat1 = zero
      src1 = zero
      srcjac = zero
      endif
	return
	end
