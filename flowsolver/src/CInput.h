/*********************************************************************

Copyright (c) 2000-2007, Stanford University, 
    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
    Charles A. Taylor (see SimVascular Acknowledgements file 
    for additional contributors to the source code).

All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions 
are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution. 
Neither the name of the Stanford University or Rensselaer Polytechnic
Institute nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior 
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

**********************************************************************/

#ifndef _CInput_H_
#define _CInput_H_
#include <vector>
#include <string>
#include <map>
#include <boost/algorithm/string/predicate.hpp>

#include "ValType.h"

using namespace std;

// http://www.cplusplus.com/reference/map/map/map/
struct caseInsensitveStringLessThanComparator {
  bool operator() (const std::string& lhs, const std::string& rhs) const
  {
    return boost::algorithm::ilexicographical_compare<std::string, std::string>(lhs, rhs);
  }
};

class CInput {
public:
  CInput(const string &, const string &default_fname = "");
  CInput(const char*, const char* = "");
  ~CInput();

  // return the entire input map
  map<string, string, caseInsensitveStringLessThanComparator> InputMap() const;

  // returns the desired string
  //  const string &GetValue(const string &) const;
  ValType GetValue(const string &) const;

  // echo the entire input map
  void EchoInputMap(const ostream &ofile);

private:

  void trim_string(string *str);

  void get_input_lines(vector<string> *, ifstream& );
  void build_map(map<string, string, caseInsensitveStringLessThanComparator> *, vector<string> *);

  map<string, string, caseInsensitveStringLessThanComparator> *input_map;
  map<string, string, caseInsensitveStringLessThanComparator> *default_map;

  vector<string> *input_text;
  vector<string> *default_text;

};




#endif
