c
c  Copyright (c) 2000-2007, Stanford University, 
c     Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
c     Charles A. Taylor (see SimVascular Acknowledgements file 
c     for additional contributors to the source code).
c
c  All rights reserved.
c
c  Redistribution and use in source and binary forms, with or without 
c  modification, are permitted provided that the following conditions 
c  are met:
c
c  Redistributions of source code must retain the above copyright notice,
c  this list of conditions and the following disclaimer. 
c  Redistributions in binary form must reproduce the above copyright 
c  notice, this list of conditions and the following disclaimer in the 
c  documentation and/or other materials provided with the distribution. 
c  Neither the name of the Stanford University or Rensselaer Polytechnic
c  Institute nor the names of its contributors may be used to endorse or
c  promote products derived from this software without specific prior 
c  written permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
c  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
c  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
c  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
c  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
c  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
c  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
c  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
c  DAMAGE.
c
c
        subroutine getFld (T,      cp,     rmu,    rlm,    rlm2mu,
     &                     con)
c
c----------------------------------------------------------------------
c
c This routine calculates the fluid material properties.
c
c input:
c  T      (npro)        : temperature
c  cp     (npro)        : specific heat at constant pressure
c
c output:
c  rmu    (npro)        : Mu
c  rlm    (npro)        : Lambda
c  rlm2mu (npro)        : Lambda + 2 Mu
c  con    (npro)        : Conductivity
c
c Note: material type flags
c         matflg(2):
c          eq. 0, constant viscosity
c          eq. 1, generalized Sutherland viscosity
c         matflg(3):
c          eq. 0, Stokes approximation
c          eq. 1, shear proportional bulk viscosity
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension T(npro),                   cp(npro),
     &            rmu(npro),                 rlm(npro),
     &            rlm2mu(npro),              con(npro)
c
c
c.... constant viscosity
c
        if (matflg(2,1) .eq. 0) then
c
          rmu = datmat(1,2,1)
c
        else
c
c.... generalized Sutherland viscosity
c
          rmu = datmat(1,2,1) * (T/datmat(2,2,1))*sqrt(T/datmat(2,2,1))
     &        * ( datmat(2,2,1) + datmat(3,2,1) ) / (T + datmat(3,2,1))
c
        endif
c
c.... calculate the second viscosity coefficient
c
        if (matflg(3,1) .eq. 0) then
c
          rlm = -pt66 * rmu
c
        else
c
          rlm = (datmat(1,3,1) - pt66) * rmu
c
        endif
c
c.... calculate the remaining quantities
c
        cp     = datmat(1,3,1)
        rlm2mu = rlm + two * rmu
        con    = datmat(1,4,1)
c
c.... return
c
        return
        end
