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

#ifndef CVSOLVERIO_H
#define CVSOLVERIO_H

#include <stdio.h>
#include <string>
#include <cstring>

#ifdef WIN32

#define openfile_ OPENFILE
#define closefile_ CLOSEFILE
#define readheader_ READHEADER
#define readdatablock_ READDATABLOCK
#define writeheader_ WRITEHEADER
#define writedatablock_ WRITEDATABLOCK
#define writestring_ WRITESTRING

#endif

#ifdef __cplusplus

class phastaIO {

public:
    
    phastaIO();
    ~phastaIO();
    
    // file control
    int openFile (const char* filename, const char *mode);
    int closeFile ();
    int rewindFile();
    
    // read info
    int readHeader (const char* keyphrase,
                  int* valueArray,
                  int   nItems,
                  const char*  datatype,
                    const char*  iotype);
    int readDataBlock (const char* keyphrase,
                       void* valueArray,
                       int  nItems,
                       const char*  datatype,
                       const char*  iotype );
    int readString(char* line);

    // write info
    int writeHeader (const char* keyphrase,
                     void* valueArray,
                     int nItems,
                     int ndataItems,
                     const char* datatype,
                     const char* iotype  );

    int writeDataBlock ( const char* keyphrase,
                         void* valueArray,
                         int nItems,
                         const char* datatype,
                         const char* iotype );
    
    int writeString (const char* string );

    // helper functions
    char* StringStripper ( const char*  istring );
    int cscompare ( const char* s1, const char* s2);
    void isBinary ( const char* iotype );
    size_t typeSize ( const char* typestring );
    void SwapArrayByteOrder ( void* array, int nbytes, int nItems );
    
private:

    FILE *filePointer_;

    std::string fileName_;
    std::string originalPreColonTokenOnLine_;
    std::string previousOriginalPreColonTokenOnLine_;

    bool byte_order_;
    int type_of_data_  ;
    int DataSize_  ;

    char LastHeaderKey_[1024];
    bool LastHeaderNotFound_;
    
    bool Wrong_Endian_;
    bool binary_format_;
    
    char *mode_;
    char *fname_;
    char Line_[1024];

};

#endif

#ifdef __cplusplus
extern "C" {
#endif
    
void openfile_( const char* filename, 
                const char* mode,
                int*  fileDescriptor );

void closefile_( int* fileDescriptor, 
                 const char* mode );

void readheader_( int* fileDescriptor,
                  const char* keyphrase,
                  void* valueArray,
                  int*  nItems,
                  const char*  datatype,
                  const char*  iotype );

void
readdatablock_( int*  fileDescriptor,
                const char* keyphrase,
                void* valueArray,
                int*  nItems,
                const char*  datatype,
                const char*  iotype );

void 
writeheader_ (  int*  fileDescriptor,
                const char* keyphrase,
                void* valueArray,
                int* nItems,
                int* ndataItems,
                const char* datatype,
                const char* iotype  );

void 
writedatablock_( int* fileDescriptor,
                 const char* keyphrase,
                 void* valueArray,
                 int* nItems,
                 const char* datatype,
                 const char* iotype );

void 
writestring_( int* fileDescriptor,
              const char* string );

void bzero_old(void* ptr, size_t sz);

#ifdef __cplusplus
}
#endif

#endif






