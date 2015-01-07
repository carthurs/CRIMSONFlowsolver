/* fortran wrapper for C function TriShapeAndDrv which returns
   the shape functions and their derivatives for tris
   */

int TriShapeAndDrv(int p,double par[2],double N[],double dN[][2]);

extern void shptri(int p, double par[], double N[], double dN[][2])
{
 TriShapeAndDrv(p,par,N,dN);

}
