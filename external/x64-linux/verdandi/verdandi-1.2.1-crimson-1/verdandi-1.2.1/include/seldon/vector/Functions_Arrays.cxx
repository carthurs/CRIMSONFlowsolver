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


#ifndef SELDON_FILE_FUNCTIONS_ARRAYS_CXX

/*
  Functions defined in this file

  QuickSort(m, n, X);
  QuickSort(m, n, X, Y);
  QuickSort(m, n, X, Y, Z);

  MergeSort(m, n, X);
  MergeSort(m, n, X, Y);
  MergeSort(m, n, X, Y, Z);

  Sort(m, n, X);
  Sort(m, n, X, Y);
  Sort(m, n, X, Y, Z);
  Sort(m, X);
  Sort(m, X, Y);
  Sort(m, X, Y, Z);
  Sort(X);
  Sort(X, Y);
  Sort(X, Y, Z);

  Assemble(m, X);
  Assemble(m, X, Y);

  RemoveDuplicate(m, X);
  RemoveDuplicate(m, X, Y);
  RemoveDuplicate(X);
  RemoveDuplicate(X, Y);
*/

namespace Seldon
{


  ////////////
  //  SORT  //


  //! Intermediary function used for quick sort algorithm.
  template<class T, class Storage, class Allocator>
  int PartitionQuickSort(int m, int n,
			 Vector<T, Storage, Allocator>& t)
  {
    T temp, v;
    v = t(m);
    int i = m - 1;
    int j = n + 1;

    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t(j) > v);

	do
	  {
	    i++;
	  }
	while (t(i) < v);

	if (i < j)
	  {
	    temp = t(i);
	    t(i) = t(j);
	    t(j) = temp;
	  }
	else
	  {
	    return j;
	  }
      }
  }


  //! Vector \a t is sorted by using QuickSort algorithm.
  /*!
    Sorts array \a t between position m and n.
  */
  template<class T, class Storage, class Allocator>
  void QuickSort(int m, int n,
		 Vector<T, Storage, Allocator>& t)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t);
	QuickSort(m, p, t);
	QuickSort(p+1, n, t);
      }
  }


  //! Intermediary function used for quick sort algorithm.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2)
  {
    T1 temp1, v;
    T2 temp2;
    v = t1(m);
    int i = m - 1;
    int j = n + 1;

    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t1(j) > v);
	do
	  {
	    i++;
	  }
	while (t1(i) < v);

	if (i < j)
	  {
	    temp1 = t1(i);
	    t1(i) = t1(j);
	    t1(j) = temp1;
	    temp2 = t2(i);
	    t2(i) = t2(j);
	    t2(j) = temp2;
	  }
	else
	  {
	    return j;
	  }
      }
  }


  //! Vector \a t1 is sorted by using QuickSort algorithm.
  /*! Sorts array \a t2 between position m and n, the sorting operation
    affects vector \a t2.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t1, t2);
	QuickSort(m, p, t1, t2);
	QuickSort(p+1, n, t1, t2);
      }
  }


  //! Intermediary function used for quick sort algorithm.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2,
			 Vector<T3, Storage3, Allocator3>& t3)
  {
    T1 temp1, v;
    T2 temp2;
    T3 temp3;
    v = t1(m);
    int i = m - 1;
    int j = n + 1;

    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t1(j) > v);

	do
	  {
	    i++;
	  }
	while (t1(i) < v);

	if (i < j)
	  {
	    temp1 = t1(i);
	    t1(i) = t1(j);
	    t1(j) = temp1;
	    temp2 = t2(i);
	    t2(i) = t2(j);
	    t2(j) = temp2;
	    temp3 = t3(i);
	    t3(i) = t3(j);
	    t3(j) = temp3;
	  }
	else
	  {
	    return j;
	  }
      }
  }


  //! Vector \a t1 is sorted by using QuickSort algorithm.
  /*! Sorts array \a t1 between position m and n, the sorting operation
    affects vectors \a t2 and \a t3.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2,
		 Vector<T3, Storage3, Allocator3>& t3)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t1, t2, t3);
	QuickSort(m, p, t1, t2, t3);
	QuickSort(p+1, n, t1, t2, t3);
      }
  }


  //! Vector \a tab1 is sorted by using MergeSort algorithm.
  /*!
    Sorts array \a tab1 between position m and n.
  */
  template<class T, class Storage, class Allocator>
  void MergeSort(int m, int n, Vector<T, Storage, Allocator>& tab1)
  {
    if (m >= n)
      return;

    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T, Storage, Allocator> tab1t(n - m + 1);
    // Performs a merge sort with a recurrence.
    // inc = 1, 2, 4, 8, ...
    while (inc < (n - m + 1) )
      {
	// 'i' is the first index of the sub array of size 2*inc.
	// A loop is performed on these sub-arrays.
	// Each sub array is divided in two sub-arrays of size 'inc'.
	// These two sub-arrays are then merged.
	for (i = m; i <= n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc; // Index of the second sub-array.
	    sup = i + 2 * inc; // End of the merged array.
	    if (sup >= n + 1)
	      sup = n + 1;

	    j = i;
	    // Loop on values of the first sub-array.
	    while (j < i + inc)
	      {
		// If the second sub-array has still unsorted elements.
		if (current < sup)
		  {
		    // Insert the elements of the second sub-array in the
		    // merged array until tab1(j) < tab1(current).
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind-m) = tab1(current);
			current++;
			ind++;
		      }

		    // Inserts the element of the first sub-array now.
		    tab1t(ind - m) = tab1(j);
		    ind++;
		    j++;

		    // If the first sub-array is sorted, all remaining
		    // elements of the second sub-array are inserted.
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind - m) = tab1(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    // If the second sub-array is sorted, all remaining
		    // elements of the first sub-array are inserted.
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind - m) = tab1(current);
			ind++;
		      }
		    j = current + 1;
		  }
	      }
	  }

	for (i = m; i < ind; i++)
	  tab1(i) = tab1t(i - m);

	inc = 2 * inc;
      }
  }


  //! Vector \a tab1 is sorted by using MergeSort algorithm.
  /*! Sorts array \a tab1 between position m and n. The sort operation affects
    \a tab2.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2)
  {
    if (m >= n)
      return;

    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T1, Storage1, Allocator1> tab1t(n - m + 1);
    Vector<T2, Storage2, Allocator2> tab2t(n - m + 1);

    while (inc < n - m + 1)
      {
	for (i = m; i <= n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc;
	    sup = i + 2 * inc;
	    if (sup >= n+1)
	      sup = n+1;

	    j = i;
	    while (j < i + inc)
	      {
		if (current < sup)
		  {
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind - m) = tab1(current);
			tab2t(ind - m) = tab2(current);
			current++;
			ind++;
		      }
		    tab1t(ind - m) = tab1(j);
		    tab2t(ind - m) = tab2(j);
		    ind++;
		    j++;
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind - m) = tab1(j);
			    tab2t(ind - m) = tab2(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind - m) = tab1(current);
			tab2t(ind - m) = tab2(current);
			ind++;
		      }
		    j = current + 1;
		  }
	      }
	  }
	for (i = m; i < ind; i++)
	  {
	    tab1(i) = tab1t(i - m);
	    tab2(i) = tab2t(i - m);
	  }
	inc = 2 * inc;
      }
  }


  //! Vector \a tab1 is sorted by using MergeSort algorithm.
  /*! Sorts array \a tab1 between position m and n. The sort operation affects
    \a tab2 and \a tab3.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2,
		 Vector<T3, Storage3, Allocator3>& tab3)
  {
    if (m >= n)
      return;

    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T1, Storage1, Allocator1> tab1t(n - m + 1);
    Vector<T2, Storage2, Allocator2> tab2t(n - m + 1);
    Vector<T3, Storage3, Allocator3> tab3t(n - m + 1);

    while (inc < n - m + 1)
      {
	for (i = m; i <= n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc;
	    sup = i + 2 * inc;
	    if (sup >= n+1)
	      sup = n+1;

	    j = i;
	    while (j < i + inc)
	      {
		if (current < sup)
		  {
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind - m) = tab1(current);
			tab2t(ind - m) = tab2(current);
			tab3t(ind - m) = tab3(current);
			current++;
			ind++;
		      }
		    tab1t(ind - m) = tab1(j);
		    tab2t(ind - m) = tab2(j);
		    tab3t(ind - m) = tab3(j);
		    ind++;
		    j++;
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind - m) = tab1(j);
			    tab2t(ind - m) = tab2(j);
			    tab3t(ind - m) = tab3(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind - m) = tab1(current);
			tab2t(ind - m) = tab2(current);
			tab3t(ind - m) = tab3(current);
			ind++;
		      }
		    j = current+1;
		  }
	      }
	  }
	for (i = m; i < ind; i++)
	  {
	    tab1(i) = tab1t(i - m);
	    tab2(i) = tab2t(i - m);
	    tab3(i) = tab3t(i - m);
	  }
	inc = 2 * inc;
      }
  }


  //! Assembles a sparse vector.
  /*! The function sorts the indices in \a Node and adds the corresponding
    values of \a Vect corresponding.

    For example, if Node = [3 2 2 0], Vect = [1.0 0.4 0.4 -0.3], the function
    will return: n = 3, Node = [0 2 3], Vect = [-0.3 0.8 1.0].

    \param[in,out] n on entry, the number of elements to assemble; on exit,
    the number of values after assembling.
    \param[in,out] Node positions in the vector.
    \param[in,out] Vect values.
    \warning Vectors \a Node and \a Vect are not resized.
  */
  template<class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2 >
  void Assemble(int& n, Vector<int, Storage1, Allocator1>& Node,
		Vector<T2, Storage2, Allocator2>& Vect)
  {
    if (n <= 1)
      return;

    Sort(n, Node, Vect);
    int prec = Node(0);
    int nb = 0;
    for (int i = 1; i < n; i++)
      if (Node(i) == prec)
	{
	  Vect(nb) += Vect(i);
	}
      else
	{
	  nb++;
	  Node(nb) = Node(i);
	  Vect(nb) = Vect(i);
	  prec = Node(nb);
	}
    n = nb + 1;
  }


  //! Sorts and removes duplicate entries of a vector.
  /*!
    \param[in,out] n on entry, the number of elements to assemble; on exit,
    the number of values after assembling.
    \param[in,out] Node vector to be assembled.
    \warning The vector \a Node is not resized.
  */
  template<class T, class Storage1, class Allocator1>
  void Assemble(int& n, Vector<T, Storage1, Allocator1>& Node)
  {
    if (n <= 1)
      return;

    Sort(n, Node);
    T prec = Node(0);
    int nb = 1;
    for (int i = 1; i < n; i++)
      if (Node(i) != prec)
	{
	  Node(nb) = Node(i);
	  prec = Node(nb);
	  nb++;
	}
    n = nb;
  }


  //! Sorts and removes duplicate entries of a vector.
  template<class T, class Storage1, class Allocator1>
  void Assemble(Vector<T, Storage1, Allocator1>& Node)
  {
    int nb = Node.GetM();
    Assemble(nb, Node);
    Node.Resize(nb);
  }


  //! Sorts and removes duplicate entries of a vector.
  /*!
    Sorting operations on \a Node also affect \a Node2.
  */
  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2)
  {
    if (n <= 1)
      return;

    Sort(n, Node, Node2);
    T prec = Node(0);
    int nb = 1;
    for (int i = 1; i < n; i++)
      if (Node(i) != prec)
	{
	  Node(nb) = Node(i);
	  Node2(nb) = Node2(i);
	  prec = Node(nb);
	  nb++;
	}
    n = nb;
  }


  //! Sorts and removes duplicate entries of a vector.
  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node)
  {
    Assemble(n, Node);
  }


  //! Sorts and removes duplicate entries of a vector.
  /*!
    Sorting operations of \a Node also affect \a Node2.
  */
  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2)
  {
    int n = Node.GetM();
    if (n <= 1)
      return;

    RemoveDuplicate(n, Node, Node2);
    Node.Resize(n);
    Node2.Resize(n);
  }


  //! Sorts and removes duplicate entries of a vector.
  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node)
  {
    int n = Node.GetM();
    if (n <= 1)
      return;

    Assemble(n, Node);
    Node.Resize(n);
  }


  //! Sorts vector \a V between a start position and an end position.
  template<class T, class Storage, class Allocator>
  void Sort(int m, int n, Vector<T, Storage, Allocator>& V)
  {
    MergeSort(m, n, V);
  }


  //! Sorts vector \a V between a start position and an end position.
  /*!
    The sorting operation of \a V also affects \a V2.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    MergeSort(m, n, V, V2);
  }


  //! Sorts vector \a V between a start position and an end position.
  /*!
    The sorting operation of \a V also affects \a V2 and \a V3.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    MergeSort(m, n, V, V2, V3);
  }


  //! Sorts the \a n first elements of \a V.
  template<class T, class Storage, class Allocator>
  void Sort(int n, Vector<T, Storage, Allocator>& V)
  {
    Sort(0, n - 1, V);
  }


  //! Sorts the \a n first elements of \a V.
  /*!
    The sorting operation of \a V also affects \a V2.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    Sort(0, n - 1, V, V2);
  }


  //! Sorts the \a n first elements of \a V.
  /*!
    The sorting operation of \a V also affects \a V2 and \a V3.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    Sort(0, n - 1, V, V2, V3);
  }


  //! Sorts vector \a V.
  template<class T, class Storage, class Allocator>
  void Sort(Vector<T, Storage, Allocator>& V)
  {
    Sort(0, V.GetM() - 1, V);
  }


  //! Sorts vector \a V.
  /*!
    The sorting operation of \a V also affects \a V2.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    Sort(0, V.GetM() - 1, V, V2);
  }


  //! Sorts vector \a V.
  /*!
    The sorting operation of \a V also affects \a V2 and \a V3.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    Sort(0, V.GetM() - 1, V, V2, V3);
  }


  //  SORT  //
  ////////////


} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_ARRAYS_CXX
#endif
