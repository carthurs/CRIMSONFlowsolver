// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_DEBUG_LEVEL_4

#include "Seldon.hxx"

using namespace Seldon;

typedef double Real_wp;
typedef complex<double> Complex_wp;


int main(int argc, char** argv)
{
  int m = 5, n = 4; IVect ipivot;
  Vector<Real_wp> xr(n), br(n), br_trans(n);
  Vector<Complex_wp> xc(n), bc(n), bc_trans(n), bc_transConj(n);
  Real_wp anorm_one, anorm_infty, rcond, ferr, berr, amax;
  Real_wp row_condition_number, col_condition_number;
  Vector<Real_wp> row_scale(m), col_scale(n);

  {
    Matrix<Real_wp, General, RowMajor> A(m,n); A.Fill(0);
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;

    DISP(A);
    A.Resize(m-1, n-1); DISP(A);
    A.Resize(m-2, n-2); DISP(A);
  }

  {
    Matrix<Real_wp, General, ColMajor> A(m,n); A.Fill(0);
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;

    DISP(A);
    A.Resize(m-1, n-1); DISP(A);
    A.Resize(m-2, n-2); DISP(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ColHerm> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;

    DISP(A); A.Resize(4,4); DISP(A);
    A.Resize(3,3); DISP(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowHerm> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;

    DISP(A); A.Resize(4,4); DISP(A);
    A.Resize(3,3); DISP(A);
  }

  /*** RowMajor ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////"<<endl;
  cout<<"//    RowMajor    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowMajor> A(n,n), Alu;
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;

    xr.Fill();
    Mlt(A, xr, br);
    MltAdd(1.0, SeldonTrans, A, xr, 0.0, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, Alu, ipivot, xr);
    RefineSolutionLU(SeldonTrans, A, Alu, ipivot, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

    Matrix<Real_wp, General, RowMajor> B(m,n);
    for (int i=0 ; i<m ; i++)
      for (int j=0 ; j<n ; j++)
	B(i,j) = i*n+j;

    GetScalingFactors(B, row_scale, col_scale, row_condition_number,
		      col_condition_number, amax);
    cout<<"Rectangular Matrix B"<<endl<<B<<endl;
    cout<<"scaling on rows "<<endl<<row_scale<<endl;
    cout<<"scaling on columns "<<endl<<col_scale<<endl;
    DISP(amax);
    DISP(row_condition_number); DISP(col_condition_number);
  }

  {
    Matrix<Complex_wp, General, RowMajor> A(n,n), Alu;
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    MltAdd(Complex_wp(1), SeldonTrans, A, xc, Complex_wp(0), bc_trans);
    MltAdd(Complex_wp(1), SeldonConjTrans, A, xc, Complex_wp(0), bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, Alu, ipivot, xc);
    RefineSolutionLU(SeldonTrans, A, Alu, ipivot, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, Alu, ipivot, xc);
    RefineSolutionLU(SeldonConjTrans, A, Alu, ipivot, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side bconj"<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = bconj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

    Matrix<Complex_wp, General, RowMajor> B(m,n);
    for (int i=0 ; i<m ; i++)
      for (int j=0 ; j<n ; j++)
	B(i,j) = i*n+j;

    GetScalingFactors(B, row_scale, col_scale, row_condition_number,
		      col_condition_number, amax);
    cout<<"Rectangular Matrix B"<<endl<<B<<endl;
    cout<<"scaling on rows "<<endl<<row_scale<<endl;
    cout<<"scaling on columns "<<endl<<col_scale<<endl;
    DISP(amax);
    DISP(row_condition_number); DISP(col_condition_number);
  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowMajor    //"<<endl;
  cout<<"////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColMajor ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////"<<endl;
  cout<<"//    ColMajor    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColMajor> A(n,n), Alu;
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;
    xr.Fill();
    Mlt(A, xr, br);
    MltAdd(1.0, SeldonTrans, A, xr, 0.0, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, Alu, ipivot, xr);
    RefineSolutionLU(SeldonTrans, A, Alu, ipivot, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

    Matrix<Real_wp, General, ColMajor> B(m,n);
    for (int i=0 ; i<m ; i++)
      for (int j=0 ; j<n ; j++)
	B(i,j) = i*n+j;

    GetScalingFactors(B, row_scale, col_scale, row_condition_number,
		      col_condition_number, amax);
    cout<<"Rectangular Matrix B"<<endl<<B<<endl;
    cout<<"scaling on rows "<<endl<<row_scale<<endl;
    cout<<"scaling on columns "<<endl<<col_scale<<endl;
    DISP(amax);
    DISP(row_condition_number); DISP(col_condition_number);
  }

  {
    Matrix<Complex_wp, General, ColMajor> A(n,n), Alu;
    A(0,0) = 2.0; A(0,1) = -1.0; A(0,2) = 0.0; A(0,3) = 1.0;
    A(1,0) = 3.0; A(1,1) = 5.0;  A(1,2) = -2.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0; A(2,3) = 0.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    MltAdd(Complex_wp(1), SeldonTrans, A, xc, Complex_wp(0), bc_trans);
    MltAdd(Complex_wp(1), SeldonConjTrans, A, xc, Complex_wp(0), bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, Alu, ipivot, xc);
    RefineSolutionLU(SeldonTrans, A, Alu, ipivot, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, Alu, ipivot, xc);
    RefineSolutionLU(SeldonConjTrans, A, Alu, ipivot, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side bconj"<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = bconj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

    Matrix<Complex_wp, General, ColMajor> B(m,n);
    for (int i=0 ; i<m ; i++)
      for (int j=0 ; j<n ; j++)
	B(i,j) = i*n+j;

    GetScalingFactors(B, row_scale, col_scale, row_condition_number,
		      col_condition_number, amax);
    cout<<"Rectangular Matrix B"<<endl<<B<<endl;
    cout<<"scaling on rows "<<endl<<row_scale<<endl;
    cout<<"scaling on columns "<<endl<<col_scale<<endl;
    DISP(amax);
    DISP(row_condition_number); DISP(col_condition_number);
  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColMajor    //"<<endl;
  cout<<"////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColSymPacked and Upper ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////"<<endl;
  cout<<"//    ColSymPacked    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, Symmetric, ColSymPacked> A(n,n), Alu;
    A(0,0) = 2.0;
    A(1,0) = 3.0; A(1,1) = 5.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;
    xr.Fill();
    Mlt(A, xr, br);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    Matrix<Complex_wp, Symmetric, ColSymPacked> A(n,n), Alu;
    A(0,0) = 2.0;
    A(1,0) = 3.0; A(1,1) = 5.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColSymPacked    //"<<endl;
  cout<<"////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowSymPacked and Upper ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////"<<endl;
  cout<<"//    RowSymPacked    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, Symmetric, RowSymPacked> A(n,n), Alu;
    A(0,0) = 2.0;
    A(1,0) = 3.0; A(1,1) = 5.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = 6.0;
    xr.Fill();
    Mlt(A, xr, br);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    Matrix<Complex_wp, Symmetric, RowSymPacked> A(n,n), Alu;
    A(0,0) = 2.0;
    A(1,0) = 3.0; A(1,1) = 5.0;
    A(2,0) = 0.0; A(2,1) = -3.0; A(2,2) = -4.0;
    A(3,0) = 2.0; A(3,1) = 4.0;  A(3,2) = 0.0;  A(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowSymPacked    //"<<endl;
  cout<<"////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColSym and Upper ***/


  cout<<endl<<endl<<endl;
  cout<<"//////////////////"<<endl;
  cout<<"//    ColSym    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, Symmetric, ColSym> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill();
    Mlt(A, xr, br);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    Matrix<Complex_wp, Symmetric, ColSym> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColSym    //"<<endl;
  cout<<"//////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowSym and Upper ***/


  cout<<endl<<endl<<endl;
  cout<<"//////////////////"<<endl;
  cout<<"//    RowSym    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, Symmetric, RowSym> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill();
    Mlt(A, xr, br);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);
    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xr);
    RefineSolutionLU(A, Alu, ipivot, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xr.Copy(br); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xr);
    cout<<"Solution of A x = b"<<endl<<xr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    Matrix<Complex_wp, Symmetric, RowSym> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = Complex_wp(6.0, 2.0);
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowSym    //"<<endl;
  cout<<"//////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColUpTriang and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<"//    ColUpTriang    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColUpTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColUpTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColUpTriang    //"<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColUpTriang and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    ColUpTriang Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColUpTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColUpTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColUpTriang Unit   //"<<endl;
  cout<<"///////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColLoTriang and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<"//    ColLoTriang    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColLoTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColLoTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColLoTriang    //"<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColLoTriang and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    ColLoTriang Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColLoTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColLoTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColLoTriang Unit   //"<<endl;
  cout<<"///////////////////////////"<<endl;
  cout<<endl<<endl<<endl;

  /*** ColUpTriangPacked and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    ColUpTriangPacked   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColUpTriangPacked> A(n,n);
    A.Get(0,0) = 2.0;
    A.Get(0,1) = 3.0; A.Get(1,1) = 5.0;
    A.Get(0,2) = 0.0; A.Get(1,2) = -3.0; A.Get(2,2) = -4.0;
    A.Get(0,3) = 2.0; A.Get(1,3) = 4.0;  A.Get(2,3) = 0.0;  A.Get(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColUpTriangPacked> A(n,n);
    A.Get(0,0) = 2.0;
    A.Get(0,1) = 3.0; A.Get(1,1) = 5.0;
    A.Get(0,2) = 0.0; A.Get(1,2) = -3.0; A.Get(2,2) = -4.0;
    A.Get(0,3) = 2.0; A.Get(1,3) = 4.0;  A.Get(2,3) = 0.0;  A.Get(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColUpTriangPacked    //"<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColUpTriangPacked and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////////"<<endl;
  cout<<"//    ColUpTriangPacked Unit  //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColUpTriangPacked> A(n,n);
    A.Get(0,0) = 1.0;
    A.Get(0,1) = 3.0; A.Get(1,1) = 1.0;
    A.Get(0,2) = 0.0; A.Get(1,2) = -3.0; A.Get(2,2) = 1.0;
    A.Get(0,3) = 2.0; A.Get(1,3) = 4.0;  A.Get(2,3) = 0.0; A.Get(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit,A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColUpTriangPacked> A(n,n);
    A.Get(0,0) = 1.0;
    A.Get(0,1) = 3.0; A.Get(1,1) = 1.0;
    A.Get(0,2) = 0.0; A.Get(1,2) = -3.0; A.Get(2,2) = 1.0;
    A.Get(0,3) = 2.0; A.Get(1,3) = 4.0;  A.Get(2,3) = Complex_wp(0.0,2.0); A.Get(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColUpTriangPacked Unit  //"<<endl;
  cout<<"////////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColLoTriangPacked and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<"//    ColLoTriangPacked    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColLoTriangPacked> A(n,n);
    A.Get(0,0) = 2.0;
    A.Get(1,0) = 3.0; A.Get(1,1) = 5.0;
    A.Get(2,0) = 0.0; A.Get(2,1) = -3.0; A.Get(2,2) = -4.0;
    A.Get(3,0) = 2.0; A.Get(3,1) = 4.0;  A.Get(3,2) = 0.0;  A.Get(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColLoTriangPacked> A(n,n);
    A.Get(0,0) = 2.0;
    A.Get(1,0) = 3.0; A.Get(1,1) = 5.0;
    A.Get(2,0) = 0.0; A.Get(2,1) = -3.0; A.Get(2,2) = -4.0;
    A.Get(3,0) = 2.0; A.Get(3,1) = 4.0;  A.Get(3,2) = 0.0;  A.Get(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColLoTriangPacked    //"<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColLoTriangPacked and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"//////////////////////////////////"<<endl;
  cout<<"//    ColLoTriangPacked Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, ColLoTriangPacked> A(n,n);
    A.Get(0,0) = 1.0;
    A.Get(1,0) = 3.0; A.Get(1,1) = 1.0;
    A.Get(2,0) = 0.0; A.Get(2,1) = -3.0; A.Get(2,2) = 1.0;
    A.Get(3,0) = 2.0; A.Get(3,1) = 4.0;  A.Get(3,2) = 0.0; A.Get(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, ColLoTriangPacked> A(n,n);
    A.Get(0,0) = 1.0;
    A.Get(1,0) = 3.0; A.Get(1,1) = 1.0;
    A.Get(2,0) = 0.0; A.Get(2,1) = -3.0; A.Get(2,2) = 1.0;
    A.Get(3,0) = 2.0; A.Get(3,1) = 4.0;  A.Get(3,2) = Complex_wp(0.0,2.0); A.Get(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColLoTriangPacked Unit  //"<<endl;
  cout<<"////////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;

  /*** RowUpTriang and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<"//    RowUpTriang    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowUpTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowUpTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowUpTriang    //"<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowUpTriang and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    RowUpTriang Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowUpTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowUpTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowUpTriang Unit   //"<<endl;
  cout<<"///////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowLoTriang and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<"//    RowLoTriang    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowLoTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowLoTriang> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowLoTriang    //"<<endl;
  cout<<"///////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowLoTriang and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    RowLoTriang Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowLoTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowLoTriang> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowLoTriang Unit    //"<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;

  /*** RowUpTriangPacked and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"////////////////////////////"<<endl;
  cout<<"//    RowUpTriangPacked   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowUpTriangPacked> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowUpTriangPacked> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowUpTriangPacked    //"<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowUpTriangPacked and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"//////////////////////////////////"<<endl;
  cout<<"//    RowUpTriangPacked Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowUpTriangPacked> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowUpTriangPacked> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 1.0;
    A.Val(0,2) = 0.0; A.Val(1,2) = -3.0; A.Val(2,2) = 1.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowUpTriangPacked Unit  //"<<endl;
  cout<<"////////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowLoTriangPacked and NonUnit ***/


  cout<<endl<<endl<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<"//    RowLoTriangPacked    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowLoTriangPacked> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = 6.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(A, br);
    Mlt(SeldonTrans, SeldonNonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xr);
    RefineSolutionLU(A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowLoTriangPacked> A(n,n);
    A.Val(0,0) = 2.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 5.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = -4.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0;  A.Val(3,3) = Complex_wp(6.0,2.0);
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(A, bc);
    Mlt(SeldonTrans, SeldonNonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonNonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(A, xc);
    RefineSolutionLU(A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonNonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowLoTriangPacked    //"<<endl;
  cout<<"/////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowLoTriangPacked and Unit ***/


  cout<<endl<<endl<<endl;
  cout<<"//////////////////////////////////"<<endl;
  cout<<"//    RowLoTriangPacked Unit    //"<<endl;
  cout<<endl<<endl<<endl;

  {
    // initialization of A
    Matrix<Real_wp, General, RowLoTriangPacked> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = 0.0; A.Val(3,3) = 1.0;
    xr.Fill(); br.Copy(xr); br_trans.Copy(xr);
    Mlt(SeldonNoTrans, SeldonUnit, A, br);
    Mlt(SeldonTrans, SeldonUnit, A, br_trans);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(br, xr);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xr, br, ferr, berr);
    cout<<"Right hand side b"<<endl<<br<<endl;
    cout<<"Solution of A x = b "<<endl<<xr<<endl;

    xr.Copy(br_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xr);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xr, br_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<br_trans<<endl;
    cout<<"Solution of A x = btrans "<<endl<<xr<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  {
    // initialization of A
    Matrix<Complex_wp, General, RowLoTriangPacked> A(n,n);
    A.Val(0,0) = 1.0;
    A.Val(1,0) = 3.0; A.Val(1,1) = 1.0;
    A.Val(2,0) = 0.0; A.Val(2,1) = -3.0; A.Val(2,2) = 1.0;
    A.Val(3,0) = 2.0; A.Val(3,1) = 4.0;  A.Val(3,2) = Complex_wp(0.0,2.0); A.Val(3,3) = 1.0;
    xc.Fill(); bc.Copy(xc); bc_trans.Copy(xc); bc_transConj.Copy(xc);
    Mlt(SeldonNoTrans, SeldonUnit, A, bc);
    Mlt(SeldonTrans, SeldonUnit, A, bc_trans);
    Mlt(SeldonConjTrans, SeldonUnit, A, bc_transConj);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNorm1);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(SeldonUnit, A, SeldonNormInf);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(SeldonNoTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    xc.Copy(bc_trans);
    SolveLU(SeldonTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonTrans, SeldonUnit, A, xc, bc_trans, ferr, berr);
    cout<<"Right hand side btrans"<<endl<<bc_trans<<endl;
    cout<<"Solution of A^t x = btrans "<<endl<<xc<<endl;

    xc.Copy(bc_transConj);
    SolveLU(SeldonConjTrans, SeldonUnit, A, xc);
    RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, xc, bc_transConj, ferr, berr);
    cout<<"Right hand side btransConj "<<endl<<bc_transConj<<endl;
    cout<<"Solution of A^h x = btransConj "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowLoTriangPacked Unit  //"<<endl;
  cout<<"////////////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColHerm ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////"<<endl;
  cout<<"//     ColHerm   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    Matrix<Complex_wp, Symmetric, ColHerm> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    DISP(A); DISP(xc);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColHerm  //"<<endl;
  cout<<"/////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** ColHermPacked ***/


  cout<<endl<<endl<<endl;
  cout<<"/////////////////////////"<<endl;
  cout<<"//     ColHermPacked   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    Matrix<Complex_wp, Symmetric, ColHermPacked> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A); DISP(anorm_infty);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    DISP(A); DISP(xc);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    ColHermPacked    //"<<endl;
  cout<<"/////////////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowHerm ***/


  cout<<endl<<endl<<endl;
  cout<<"///////////////////"<<endl;
  cout<<"//     RowHerm   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    Matrix<Complex_wp, Symmetric, RowHerm> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    DISP(A); DISP(xc);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//    RowHerm  //"<<endl;
  cout<<"/////////////////"<<endl;
  cout<<endl<<endl<<endl;


  /*** RowHermPacked ***/


  cout<<endl<<endl<<endl;
  cout<<"/////////////////////////"<<endl;
  cout<<"//     RowHermPacked   //"<<endl;
  cout<<endl<<endl<<endl;

  {
    Matrix<Complex_wp, Symmetric, RowHermPacked> A(n,n), Alu;
    A.Val(0,0) = 2.0;
    A.Val(0,1) = 3.0; A.Val(1,1) = 5.0;
    A.Val(0,2) = Complex_wp(-1.0,1.0); A.Val(1,2) = -3.0; A.Val(2,2) = -4.0;
    A.Val(0,3) = 2.0; A.Val(1,3) = 4.0;  A.Val(2,3) = Complex_wp(0.0,2.0);  A.Val(3,3) = 6.0;
    xc.Fill();

    Mlt(A, xc, bc);
    cout<<"Matrix A"<<endl<<A<<endl;
    Copy(bc, xc);
    anorm_one = Norm1(A);
    cout<<"1-norm of A "<<anorm_one<<endl;
    anorm_infty = NormInf(A); DISP(anorm_infty);
    cout<<"infinity-norm of A "<<anorm_infty<<endl;

    Alu.Copy(A);

    GetLU(Alu, ipivot);
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNorm1, anorm_one);
    cout<<"The reciprocal of condition number in 1-norm is equal to "<<rcond<<endl;
    rcond = ReciprocalConditionNumber(Alu, ipivot, SeldonNormInf, anorm_infty);
    cout<<"The reciprocal of condition number in infinity-norm is equal to "<<rcond<<endl;

    SolveLU(Alu, ipivot, xc);
    RefineSolutionLU(A, Alu, ipivot, xc, bc, ferr, berr);
    cout<<"Right hand side b"<<endl<<bc<<endl;
    cout<<"Solution of A x = b "<<endl<<xc<<endl;

    cout<<"Forward error "<<ferr<<endl; cout<<"Backward error "<<berr<<endl;

    xc.Copy(bc); Alu.Copy(A);
    DISP(A); DISP(xc);
    GetAndSolveLU(Alu, ipivot, xc);
    cout<<"Solution of A x = b"<<endl<<xc<<endl;

    GetInverse(A);
    cout<<"Inverse of A is equal to "<<endl<<A<<endl;

  }

  cout<<endl<<endl<<endl;
  cout<<"//     RowHermPacked   //"<<endl;
  cout<<"/////////////////////////"<<endl;
  cout<<endl<<endl<<endl;

  return 0;
}
