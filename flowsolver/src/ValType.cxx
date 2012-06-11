#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "ValType.h"

ValType::operator int()
{
  used = true;
  return get_int(str);
}

ValType::operator vector<double>()
{
  used = true;
  return get_vector(str);
}

ValType::operator vector<int>()
{
  used = true;
  return get_ivector(str);
}

ValType::operator double()
{
  used = true;
  return get_double(str);
}

ValType::operator double*()
{
  used = true;
  return get_double_array(str);
}

ValType::operator string()
{
  used = true;
  return get_string(str);
}

//
//  function implementations for type specific conversions
//

int ValType::get_int(string str)
{
  int i = atoi(str.c_str());
  return i;
}

double ValType::get_double(string str)
{
  double x = atof(str.c_str());
  return x;
}

double *ValType::get_double_array(string str)
{
  std::istringstream ist(str);
  vector<double> vec;
  double v;
  while ( ist >> v ) {
    vec.push_back(v);
  }
  int n = vec.size();
  double *x = new double[n];
  for (int i=0; i < n; i++) {
    x[i] = vec[i];
  }
  return x;
}

vector<double> ValType::get_vector(string str)
{
  std::istringstream ist(str);
  vector<double> vec;
  double v;
  while ( ist >> v ) {
    vec.push_back(v);
  }
  return vec;
}

vector<int> ValType::get_ivector(string str)
{
  std::istringstream ist(str);
  vector<int> vec;
  int v;
  while ( ist >> v ) {
    vec.push_back(v);
  }
  return vec;
}

string ValType::get_string(string str)
{
  return str;
}



