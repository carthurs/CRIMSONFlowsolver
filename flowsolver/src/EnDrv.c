/*********************************************************************

Copyright (c) 2000-2007, Stanford University, 
    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
    Charles A. Taylor (see SimVascular Acknowledgements file 
    for additional contributors to the source code).

All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions 
are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution. 
Neither the name of the Stanford University or Rensselaer Polytechnic
Institute nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior 
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

**********************************************************************/

/*

This function calculates the derivative of edge mode for simplex.

*/

#ifdef __cplusplus
extern "C" {
#endif

int EnDrv(int ip, double r, double s, double drv[2]) {
   double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
   double t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35;

   /* p=2 */
   if( ip==0 ) {
      drv[0] = 0.0;
      drv[1] = 0.0;
   /* p=3 */
   } else if( ip==1 ) {
      drv[0] = -1.0;
      drv[1] = 1.0;
   /* p=4 */
   } else if( ip==2 ) {
      drv[0] = -3.0*s+2.0*r;
      drv[1] = 2.0*s-3.0*r;
   /* p=5 */
   } else if( ip==3 ) {
      t1 = s*s;
      t2 = r*s;
      t3 = r*r;
      drv[0] = -6.0*t1+12.0*t2-3.0*t3;
      drv[1] = 3.0*t1-12.0*t2+6.0*t3;
   /* p=6 */
   } else if( ip==4 ) {
      t1 = s*s;
      t2 = t1*s;
      t3 = r*t1;
      t4 = r*r;
      t5 = t4*s;
      t6 = t4*r;
      drv[0] = -10.0*t2+40.0*t3-30.0*t5+4.0*t6;
      drv[1] = 4.0*t2-30.0*t3+40.0*t5-10.0*t6;
   /* p=7 */
   } else if( ip==5 ) {
      t1 = s*s;
      t2 = t1*t1;
      t4 = r*t1*s;
      t5 = r*r;
      t6 = t5*t1;
      t8 = t5*r*s;
      t9 = t5*t5;
      drv[0] = -15.0*t2+100.0*t4-150.0*t6+60.0*t8-5.0*t9;
      drv[1] = 5.0*t2-60.0*t4+150.0*t6-100.0*t8+15.0*t9;
   /* p=8 */
   } else if( ip==6 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*s;
      t4 = r*t2;
      t5 = r*r;
      t7 = t5*t1*s;
      t9 = t5*r*t1;
      t10 = t5*t5;
      t11 = t10*s;
      t12 = t10*r;
      drv[0] = -21.0*t3+210.0*t4-525.0*t7+420.0*t9-105.0*t11+6.0*t12;
      drv[1] = 6.0*t3-105.0*t4+420.0*t7-525.0*t9+210.0*t11-21.0*t12;
   /* p=9 */
   } else if( ip==7 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t1;
      t5 = r*t2*s;
      t6 = r*r;
      t7 = t6*t2;
      t10 = t6*r*t1*s;
      t11 = t6*t6;
      t12 = t11*t1;
      t14 = t11*r*s;
      t15 = t11*t6;
      drv[0] = -28.0*t3+392.0*t5-1470.0*t7+1960.0*t10-980.0*t12+168.0*t14-7.0*
t15;
      drv[1] = 7.0*t3-168.0*t5+980.0*t7-1960.0*t10+1470.0*t12-392.0*t14+28.0*
t15;
   /* p=10 */
   } else if( ip==8 ) {
      t1 = s*s;
      t2 = t1*s;
      t3 = t1*t1;
      t4 = t3*t2;
      t6 = r*t3*t1;
      t7 = r*r;
      t9 = t7*t3*s;
      t10 = t7*r;
      t11 = t10*t3;
      t12 = t7*t7;
      t13 = t12*t2;
      t15 = t12*r*t1;
      t17 = t12*t7*s;
      t18 = t12*t10;
      drv[0] = -36.0*t4+672.0*t6-3528.0*t9+7056.0*t11-5880.0*t13+2016.0*t15
-252.0*t17+8.0*t18;
      drv[1] = 8.0*t4-252.0*t6+2016.0*t9-5880.0*t11+7056.0*t13-3528.0*t15+672.0
*t17-36.0*t18;
   /* p=11 */
   } else if( ip==9 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t1*s;
      t6 = r*t2*t4;
      t7 = r*r;
      t9 = t7*t2*t1;
      t10 = t7*r;
      t12 = t10*t2*s;
      t13 = t7*t7;
      t14 = t13*t2;
      t16 = t13*r*t4;
      t18 = t13*t7*t1;
      t20 = t13*t10*s;
      t21 = t13*t13;
      drv[0] = -45.0*t3+1080.0*t6-7560.0*t9+21168.0*t12-26460.0*t14+15120.0*t16
-3780.0*t18+360.0*t20-9.0*t21;
      drv[1] = 9.0*t3-360.0*t6+3780.0*t9-15120.0*t12+26460.0*t14-21168.0*t16+
7560.0*t18-1080.0*t20+45.0*t21;
   /* p=12 */
   } else if( ip==10 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t3*s;
      t5 = r*t3;
      t6 = r*r;
      t7 = t1*s;
      t9 = t6*t2*t7;
      t10 = t6*r;
      t12 = t10*t2*t1;
      t13 = t6*t6;
      t15 = t13*t2*s;
      t17 = t13*r*t2;
      t19 = t13*t6*t7;
      t21 = t13*t10*t1;
      t22 = t13*t13;
      t23 = t22*s;
      t24 = t22*r;
      drv[0] = -55.0*t4+1650.0*t5-14850.0*t9+55440.0*t12-97020.0*t15+83160.0*
t17-34650.0*t19+6600.0*t21-495.0*t23+10.0*t24;
      drv[1] = 10.0*t4-495.0*t5+6600.0*t9-34650.0*t12+83160.0*t15-97020.0*t17+
55440.0*t19-14850.0*t21+1650.0*t23-55.0*t24;
   /* p=13 */
   } else if( ip==11 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t3*t1;
      t6 = r*t3*s;
      t7 = r*r;
      t8 = t7*t3;
      t9 = t7*r;
      t10 = t1*s;
      t12 = t9*t2*t10;
      t13 = t7*t7;
      t15 = t13*t2*t1;
      t18 = t13*r*t2*s;
      t20 = t13*t7*t2;
      t22 = t13*t9*t10;
      t23 = t13*t13;
      t24 = t23*t1;
      t26 = t23*r*s;
      t27 = t23*t7;
      t28 = -66.0*t4+2420.0*t6-27225.0*t8+130680.0*t12-304920.0*t15+365904.0*
t18-228690.0*t20+72600.0*t22-10890.0*t24+660.0*t26-11.0*t27;
      t29 = 11.0*t4-660.0*t6+10890.0*t8-72600.0*t12+228690.0*t15-365904.0*t18+
304920.0*t20-130680.0*t22+27225.0*t24-2420.0*t26+66.0*t27;
      drv[0] = t28;
      drv[1] = t29;
   /* p=14 */
   } else if( ip==12 ) {
      t1 = s*s;
      t2 = t1*s;
      t3 = t1*t1;
      t4 = t3*t3;
      t5 = t4*t2;
      t7 = r*t4*t1;
      t8 = r*r;
      t10 = t8*t4*s;
      t11 = t8*r;
      t12 = t11*t4;
      t13 = t8*t8;
      t15 = t13*t3*t2;
      t18 = t13*r*t3*t1;
      t21 = t13*t8*t3*s;
      t23 = t13*t11*t3;
      t24 = t13*t13;
      t25 = t24*t2;
      t27 = t24*r*t1;
      t29 = t24*t8*s;
      t30 = t24*t11;
      t31 = -78.0*t5+3432.0*t7-47190.0*t10+283140.0*t12-849420.0*t15+1359072.0*
t18-1189188.0*t21+566280.0*t23-141570.0*t25+17160.0*t27-858.0*t29+12.0*t30;
      t32 = 12.0*t5-858.0*t7+17160.0*t10-141570.0*t12+566280.0*t15-1189188.0*
t18+1359072.0*t21-849420.0*t23+283140.0*t25-47190.0*t27+3432.0*t29-78.0*t30;
      drv[0] = t31;
      drv[1] = t32;
   /* p=15 */
   } else if( ip==13 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t3*t2;
      t5 = t1*s;
      t7 = r*t3*t5;
      t8 = r*r;
      t10 = t8*t3*t1;
      t11 = t8*r;
      t13 = t11*t3*s;
      t14 = t8*t8;
      t15 = t14*t3;
      t18 = t14*r*t2*t5;
      t21 = t14*t8*t2*t1;
      t24 = t14*t11*t2*s;
      t25 = t14*t14;
      t26 = t25*t2;
      t28 = t25*r*t5;
      t30 = t25*t8*t1;
      t32 = t25*t11*s;
      t33 = t25*t14;
      t34 = -91.0*t4+4732.0*t7-78078.0*t10+572572.0*t13-2147145.0*t15+4416984.0
*t18-5153148.0*t21+3435432.0*t24-1288287.0*t26+260260.0*t28-26026.0*t30+1092.0*
t32-13.0*t33;
      t35 = 13.0*t4-1092.0*t7+26026.0*t10-260260.0*t13+1288287.0*t15-3435432.0*
t18+5153148.0*t21-4416984.0*t24+2147145.0*t26-572572.0*t28+78078.0*t30-4732.0*
t32+91.0*t33;
      drv[0] = t34;
      drv[1] = t35;
    } else
      return 0;
    return 1;
}
   
#ifdef __cplusplus
}
#endif
