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


#ifndef SELDON_FILE_COMMON_CXX

#include "Common.hxx"

template <class T>
void PrintArray(T* v, int lgth)
{
  for (int k = 0; k < lgth - 1; k++)
    std::cout << v[k] << " | ";
  std::cout << v[lgth - 1] << std::endl;
}

namespace Seldon
{


  //! Converts most types to string.
  /*!
    \param input variable to be converted.
    \return A string containing 'input'.
  */
  template<typename T>
  std::string to_str(const T& input)
  {
    std::ostringstream output;
    output << input;
    return output.str();
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \param[out] num \a s converted to 'T'.
  */
  template <class T>
  void to_num(std::string s, T& num)
  {
    std::istringstream str(s);
    str >> num;
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \return \a s converted to 'T'.
  */
  template <class T>
  T to_num(std::string s)
  {
    T num;
    std::istringstream str(s);
    str >> num;
    return num;
  }


  //! Sets a number to zero.
  /*!
    \param[in,out] number number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(T& number)
  {
    number = T(0);
  }


  //! Sets a complex number to zero.
  /*!
    \param[in,out] number complex number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(complex<T>& number)
  {
    number = complex<T>(T(0), T(0));
  }


  //! Sets a number to one.
  /*!
    \param[in,out] number number to be set to one.
  */
  template <class T>
  inline void SetComplexOne(T& number)
  {
    number = T(1);
  }


  //! Sets a complex number to (1, 0).
  /*!
    \param[in,out] number complex number to be set to (1, 0).
  */
  template <class T>
  inline void SetComplexOne(complex<T>& number)
  {
    number = complex<T>(T(1), T(0));
  }


}  // namespace Seldon.

#define SELDON_FILE_COMMON_CXX
#endif
