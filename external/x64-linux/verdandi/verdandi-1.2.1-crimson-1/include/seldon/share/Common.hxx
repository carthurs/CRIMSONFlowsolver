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


#ifndef SELDON_FILE_COMMON_HXX

#include <complex>
#include <iostream>

template <class T>
void PrintArray(T* v, int lgth);

namespace Seldon
{


  using namespace std;


  template<typename T>
  std::string to_str(const T& input);

  template <class T>
  void to_num(std::string s, T& num);

  template <class T>
  T to_num(std::string s);

  template <class T>
  void SetComplexZero(T& number);

  template <class T>
  void SetComplexZero(complex<T>& number);

  template <class T>
  void SetComplexOne(T& number);

  template <class T>
  void SetComplexOne(complex<T>& number);


}  // namespace Seldon.

#define SELDON_FILE_COMMON_HXX
#endif
