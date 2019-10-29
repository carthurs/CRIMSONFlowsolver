// Copyright (C) 2001-2010 Vivien Mallet
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


#ifndef SELDON_FILE_SHARE_MATRIXFLAG_HXX


namespace Seldon
{


  /////////////////////
  // SELDONTRANSPOSE //
  /////////////////////


  class SeldonTranspose
  {
#ifdef SELDON_WITH_BLAS
  protected:
    CBLAS_TRANSPOSE cblas_status_;
#endif

  protected:
    // 0: Trans, 1: NoTrans, 2: ConjTrans.
    int status_;

  public:
    SeldonTranspose(int status);
#ifdef SELDON_WITH_BLAS
    SeldonTranspose(const enum CBLAS_TRANSPOSE status);
    operator CBLAS_TRANSPOSE() const;
#endif
    char Char() const;
    char RevChar() const;
    bool Trans() const;
    bool NoTrans() const;
    bool ConjTrans() const;
  };


  class class_SeldonTrans: public SeldonTranspose
  {
  public:
    class_SeldonTrans();
  };


  class class_SeldonNoTrans: public SeldonTranspose
  {
  public:
    class_SeldonNoTrans();
  };


  class class_SeldonConjTrans: public SeldonTranspose
  {
  public:
    class_SeldonConjTrans();
  };


  extern class_SeldonTrans SeldonTrans;
  extern class_SeldonNoTrans SeldonNoTrans;
  extern class_SeldonConjTrans SeldonConjTrans;


  ////////////////
  // SELDONDIAG //
  ////////////////


  class SeldonDiag
  {
#ifdef SELDON_WITH_BLAS
  protected:
    CBLAS_DIAG cblas_status_;
#endif

  protected:
    // 0: NonUnit, 1: Unit.
    int status_;

  public:
    SeldonDiag(int status);
#ifdef SELDON_WITH_BLAS
    operator CBLAS_DIAG() const;
#endif
    char Char() const;
    bool NonUnit() const;
    bool Unit() const;
  };


  class class_SeldonNonUnit: public SeldonDiag
  {
  public:
    class_SeldonNonUnit();
  };


  class class_SeldonUnit: public SeldonDiag
  {
  public:
    class_SeldonUnit();
  };


  extern class_SeldonNonUnit SeldonNonUnit;
  extern class_SeldonUnit SeldonUnit;


  ////////////////
  // SELDONUPLO //
  ////////////////


  class SeldonUplo
  {
#ifdef SELDON_WITH_BLAS
  protected:
    CBLAS_UPLO cblas_status_;
#endif

  protected:
    // 0: Upper, 1: Lower.
    int status_;

  public:
    SeldonUplo(int status);
#ifdef SELDON_WITH_BLAS
    operator CBLAS_UPLO() const;
#endif
    operator char() const;
    bool Upper() const;
    bool Lower() const;
    char Char() const;
    char RevChar() const;
  };


  extern SeldonUplo SeldonUpper;
  extern SeldonUplo SeldonLower;


  ////////////////
  // SELDONNORM //
  ////////////////


  class SeldonNorm
  {
  protected:
    // 0: Infinity-norm, 1: 1-norm.
    int status_;

  public:
    SeldonNorm(int status);
    operator char() const;
    char Char() const;
    char RevChar() const;
  };


  extern SeldonNorm SeldonNormInf;
  extern SeldonNorm SeldonNorm1;


  /////////////////////
  // SELDONCONJUGATE //
  /////////////////////


  class SeldonConjugate
  {
  protected:
    // false: Unconj, true: Conj.
    bool status_;

  public:
    SeldonConjugate(bool status);
    inline bool Conj() const;
  };


  extern SeldonConjugate SeldonUnconj;
  extern SeldonConjugate SeldonConj;


  ////////////////
  // SELDONSIDE //
  ////////////////


  class SeldonSide
  {
#ifdef SELDON_WITH_BLAS
  protected:
    CBLAS_SIDE cblas_status_;
#endif

  protected:
    // 0: Left, 1: Right.
    int status_;

  public:
    SeldonSide(int status);
#ifdef SELDON_WITH_BLAS
    SeldonSide(const enum CBLAS_SIDE status);
    operator CBLAS_SIDE() const;
#endif
    bool Left() const;
    bool Right() const;
  };


  class class_SeldonLeft: public SeldonSide
  {
  public:
    class_SeldonLeft();
  };


  class class_SeldonRight: public SeldonSide
  {
  public:
    class_SeldonRight();
  };


  extern class_SeldonLeft SeldonLeft;
  extern class_SeldonRight SeldonRight;


}


#define SELDON_FILE_SHARE_MATRIXFLAG_HXX
#endif
