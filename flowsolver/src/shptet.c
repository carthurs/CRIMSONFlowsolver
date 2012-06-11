/* fortran wrapper for C function TetShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */

int TetShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

void shptet(int p, double par[], double N[], double dN[][3])
{
   TetShapeAndDrv(p,par,N,dN);
}
