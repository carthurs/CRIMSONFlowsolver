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

#include <fstream>
#include "Input.h"
#include "ValType.h"
#include <stdexcept>
#include <sstream>
#include <iostream>

// return a given key value (if it's in the map)
ValType Input::GetValue(const string &str) const
{
  if (input_map->find(str) != input_map->end()) {
    if ( (*input_map)[str] == "NODEFAULT" ) {
      std::ostringstream ost;
      ost << "required input variable not set: " << str << ends;
      throw invalid_argument( ost.str() );
    }
  } else {
    std::ostringstream ost;
    ost << "required input variable not set: " << str << ends;
    throw invalid_argument( ost.str() );
  }

  return ValType( (*input_map)[str] );
}


Input::Input(const string &fname, const string &default_fname)
{
  // open the input file
  ifstream infile( fname.c_str(), ios::in);

  if(!infile || infile.eof()){
    cerr<<" Input file does not exist or is empty or perhaps you forgot mpirun? "<<endl;
    exit(-2);
  }

  // allocate memory
  input_text   = new vector<string>;
  input_map    = new map<string,string>;

  // get the lines of text from the input file
  get_input_lines(input_text, infile);
  build_map(input_map, input_text);

  // build and merge with default map ( if desired )
  if (!default_fname.empty()) {
    ifstream infile2( default_fname.c_str(), ios::in);

    map<string,string> *default_map  = new map<string,string>;
    vector<string> *default_text = new vector<string>;

    get_input_lines(default_text, infile2);
    build_map(default_map, default_text);

    // merge the two maps
    map<string,string>::const_iterator iter = default_map->begin();
    for ( ; iter != default_map->end(); ++iter ) {
      string defkey = iter->first;
      string defval = iter->second;
      if ( input_map->find(defkey) == input_map->end() ) {
	(*input_map)[defkey] = defval;
      }
    }
    infile2.close();

    delete default_map;
    delete default_text;
    
  } else  {
    cerr << "Input warning: no input.config file found." << endl;
    cerr << "Get one from source directory." << endl;
    exit(-2);
  }

  infile.close();
  
}
  
Input::~Input()
{
  delete input_text;
  delete input_map;
}


// return the input map
map<string,string> Input::InputMap() const
{
  return *input_map;
}

// echo the entire map
void Input::EchoInputMap(const ostream &ofile)
{
  map<string,string>::const_iterator iter = input_map->begin();
  for ( ; iter != input_map->end(); ++iter ) {
    cout << "Keyphrase: " << iter->first << endl 
	 << "Keyvalue:  " << iter->second << endl << endl;
  }
}

// read the input text from the given stream
void Input::get_input_lines(vector<string> *text, ifstream &infile)
{
  string textline;
  while ( getline( infile, textline, '\n' ) ) {
    // ignore everything on a comment line
    if ( textline[0] != '#' ) {
      text->push_back( textline );
    }
  }
}


// 
void Input::build_map(map<string,string> *inmap,
		      vector<string>     *intext)
{
  // iterate through input_text of text and separate at :'s
  for (int i = 0 ; i < intext->size(); i++) {
    string textlineALL = (*intext)[i];
    string textline;
    int pos  = 0;
     
    // modification introduced so that comments starting midway in a file 
    // can be handled.

    if ( (pos = textlineALL.find_first_of( '#',pos)) != string::npos) {
      textline = textlineALL.substr(0,pos);
    }else {
      textline = textlineALL;
    }
    pos = 0;
    if ( (pos = textline.find_first_of( ':',pos)) != string::npos) {
      
      // get the keyphrase
      string keywd = textline.substr(0,pos);
      trim_string(&keywd);
      
      // get the key-value
      string keyval = textline.substr( pos+1, textline.length() - pos);
      trim_string(&keyval);
      
      // put the pair into the map
      (*inmap)[keywd] = keyval;
      
    }
  }
}

// remove leading and trailing spaces (or tabs)
void Input::trim_string(string *str)
{
  // check for empty string
  int length = str->length();
  if ( length == 0 )
    return;
  
  // erase leading spaces (or tabs)
  int pos0 = 0;
  while ( (*str)[pos0] == ' ' || (*str)[pos0] == '\t') {
    pos0++;
  }
  if ( pos0 > 0 ) {
    str->erase(0,pos0);
  }

  length = str->length();
  if ( length == 0 )
    return;

  pos0 = length-1;
  // erase trailing spaces (or tabs)
  while ( (*str)[pos0] == ' ' || (*str)[pos0] == '\t') {
    pos0--;
  }
  if ( pos0 < length-1 ) {
    str->erase(pos0+1, length-pos0);
  }

}
