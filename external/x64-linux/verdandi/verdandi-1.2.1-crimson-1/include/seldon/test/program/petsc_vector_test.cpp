// Copyright (C) 2010, INRIA
// Author(s): Marc Fragu
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



static char help[] = "Tests solving linear system on 0 by 0 matrix.\n\n";

#include <petscmat.h>
#include <petscvec.h>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

#include "vector/PetscVector.cxx"

#include <stdlib.h>
#include <petscmat.h>
#include <petscksp.h>


int main(int argc, char** args)
{

  int ierr;

  PetscInitialize(&argc, &args, (char *)0, help);

  {
    Vec x, y, z;

    int N = 10;

    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &z); CHKERRQ(ierr);
    ierr = VecSet(x, 3.0); CHKERRQ(ierr);
    ierr = VecSet(y, 1.0); CHKERRQ(ierr);
    ierr = VecSet(z, 2.0); CHKERRQ(ierr);

    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(z); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(z); CHKERRQ(ierr);

    ierr = VecAXPY(y, -1.0, x); CHKERRQ(ierr);

    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);
    ierr = VecDestroy(&z); CHKERRQ(ierr);
  }

  // Simple example with sequential PETSc Vector.
  {
    Vector<double, PETScSeq> x(10), y(10);
    x.Reallocate(10);
    y.Reallocate(10);

    x.Fill(3.);
    y.Fill(1.);

    x.Print();
    y.Print();

    Add(-1., x, y);

    y.Print();

    Mlt(2., y);

    y.Print();
  }

  // The same example with parallel PETSc Vector.
  {
    Vector<double, PETScPar> x, y;
    x.Reallocate(10);
    y.Reallocate(10);
    x.Fill(3.);
    y.Fill(1.);

    x.Print();
    y.Print();

    Add(-1., x, y);

    y.Print();

    Mlt(2., y);

    y.Print();

    y *= 2;

    y.Print();

    y.Fill();

    y.Print();
  }

  // Combined use of PETSc native variable and PETSc Vector from Seldon.
  {
    Vec z;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 10, &z); CHKERRQ(ierr);
    Vector<double, PETScSeq> x, y(z);
    x.Reallocate(10);
    x.Fill(3.);
    y.Fill(2.0);

    x.Print();

    Add(-1., y, x);

    x.Print();

    x *= 2;

    x.Print();

    x.Fill();

    x.Print();

    ierr = VecDestroy(&z); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;

}
