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


#ifndef SELDON_FILE_SHARE_MATRIXFLAG_CXX


#include "MatrixFlag.hxx"


namespace Seldon
{


  /////////////////////
  // SELDONTRANSPOSE //
  /////////////////////


  SeldonTranspose::SeldonTranspose(int status)
  {
    status_ = status;
#ifdef SELDON_WITH_BLAS
    if (status_ == 0)
      cblas_status_ = CblasTrans;
    else if (status_ == 1)
      cblas_status_ = CblasNoTrans;
    else
      cblas_status_ = CblasConjTrans;
#endif
  }


#ifdef SELDON_WITH_BLAS
  SeldonTranspose::SeldonTranspose(const enum CBLAS_TRANSPOSE status):
    cblas_status_(status)
  {
    if (cblas_status_ == CblasTrans)
      status_ = 0;
    else if (cblas_status_ == CblasNoTrans)
      status_ = 1;
    else
      status_ = 2;
  }
#endif


#ifdef SELDON_WITH_BLAS
  SeldonTranspose::operator CBLAS_TRANSPOSE() const
  {
    return cblas_status_;
  }
#endif


  char SeldonTranspose::Char() const
  {
    if (status_ == 0)
      return 'T';
    else if (status_ == 1)
      return 'N';
    else
      return 'C';
  }


  char SeldonTranspose::RevChar() const
  {
    if (status_ == 0)
      return 'N';
    else if (status_ == 1)
      return 'T';
    else
      return 'N';
  }


  bool SeldonTranspose::Trans() const
  {
    return (status_ == 0);
  }


  bool SeldonTranspose::NoTrans() const
  {
    return (status_ == 1);
  }


  bool SeldonTranspose::ConjTrans() const
  {
    return (status_ == 2);
  }


  class_SeldonTrans::class_SeldonTrans(): SeldonTranspose(0)
  {
  }


  class_SeldonNoTrans::class_SeldonNoTrans(): SeldonTranspose(1)
  {
  }


  class_SeldonConjTrans::class_SeldonConjTrans(): SeldonTranspose(2)
  {
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  class_SeldonTrans SeldonTrans;
  class_SeldonNoTrans SeldonNoTrans;
  class_SeldonConjTrans SeldonConjTrans;
#endif


  ////////////////
  // SELDONDIAG //
  ////////////////


  SeldonDiag::SeldonDiag(int status)
  {
    status_ = status;
#ifdef SELDON_WITH_BLAS
    if (status_ == 0)
      cblas_status_ = CblasNonUnit;
    else
      cblas_status_ = CblasUnit;
#endif
  }


#ifdef SELDON_WITH_BLAS
  SeldonDiag::operator CBLAS_DIAG() const
  {
    return cblas_status_;
  }
#endif


  char SeldonDiag::Char() const
  {
    return (status_ == 0) ? 'N' : 'U';
  }


  bool SeldonDiag::NonUnit() const
  {
    return (status_ == 0);
  }


  bool SeldonDiag::Unit() const
  {
    return (status_ == 1);
  }


  class_SeldonNonUnit::class_SeldonNonUnit(): SeldonDiag(0)
  {
  }


  class_SeldonUnit::class_SeldonUnit(): SeldonDiag(1)
  {
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  class_SeldonNonUnit SeldonNonUnit;
  class_SeldonUnit SeldonUnit;
#endif


  ////////////////
  // SELDONUPLO //
  ////////////////


  SeldonUplo::SeldonUplo(int status)
  {
    status_ = status;
#ifdef SELDON_WITH_BLAS
    if (status_ == 0)
      cblas_status_ = CblasUpper;
    else
      cblas_status_ = CblasLower;
#endif
  }


#ifdef SELDON_WITH_BLAS
  SeldonUplo::operator CBLAS_UPLO() const
  {
    return cblas_status_;
  }
#endif


  SeldonUplo::operator char() const
  {
    return (status_ == 0) ? 'U' : 'L';
  }


  bool SeldonUplo::Upper() const
  {
    return (status_ == 0);
  }


  bool SeldonUplo::Lower() const
  {
    return (status_ == 1);
  }


  char SeldonUplo::Char() const
  {
    return (status_ == 0) ? 'U' : 'L';
  }


  char SeldonUplo::RevChar() const
  {
    return (status_ == 0) ? 'L' : 'U';
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  SeldonUplo SeldonUpper(0);
  SeldonUplo SeldonLower(1);
#endif


  ////////////////
  // SELDONNORM //
  ////////////////


  SeldonNorm::SeldonNorm(int status)
  {
    status_ = status;
  }


  SeldonNorm::operator char() const
  {
    return (status_ == 0) ? 'I' : '1';
  }


  char SeldonNorm::Char() const
  {
    return (status_ == 0) ? 'I' : '1';
  }


  char SeldonNorm::RevChar() const
  {
    return (status_ == 0) ? '1' : 'I';
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  SeldonNorm SeldonNormInf(0);
  SeldonNorm SeldonNorm1(1);
#endif


  /////////////////////
  // SELDONCONJUGATE //
  /////////////////////


  SeldonConjugate::SeldonConjugate(bool status)
  {
    status_ = status;
  }


  inline bool SeldonConjugate::Conj() const
  {
    return status_;
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  SeldonConjugate SeldonUnconj(false);
  SeldonConjugate SeldonConj(true);
#endif


  ////////////////
  // SELDONSIDE //
  ////////////////


  SeldonSide::SeldonSide(int status)
  {
    status_ = status;
#ifdef SELDON_WITH_BLAS
    if (status_ == 0)
      cblas_status_ = CblasLeft;
    else
      cblas_status_ = CblasRight;
#endif
  }


#ifdef SELDON_WITH_BLAS
  SeldonSide::SeldonSide(const enum CBLAS_SIDE status):
    cblas_status_(status)
  {
    if (cblas_status_ == CblasLeft)
      status_ = 0;
    else
      status_ = 1;
  }
#endif


#ifdef SELDON_WITH_BLAS
  SeldonSide::operator CBLAS_SIDE() const
  {
    return cblas_status_;
  }
#endif


  bool SeldonSide::Left() const
  {
    return (status_ == 0);
  }


  bool SeldonSide::Right() const
  {
    return (status_ == 1);
  }


  class_SeldonLeft::class_SeldonLeft(): SeldonSide(0)
  {
  }


  class_SeldonRight::class_SeldonRight(): SeldonSide(1)
  {
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  class_SeldonLeft SeldonLeft;
  class_SeldonRight SeldonRight;
#endif


}


#define SELDON_FILE_SHARE_MATRIXFLAG_CXX
#endif
