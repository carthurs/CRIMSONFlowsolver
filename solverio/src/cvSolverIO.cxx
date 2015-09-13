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
#include "cvSolverIO.h"

#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <ctype.h>


#define INT 1
#define FLOAT 2
#define DOUBLE 3
#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

#define PHASTA_MAGIC_NUMBER 362436

phastaIO** phastaIOfp = NULL;

phastaIO::phastaIO () {
    //byte_order_;
    type_of_data_=0;
    DataSize_=0;
    LastHeaderNotFound_ = false;
    Wrong_Endian_ = false ;
    binary_format_ = true; 
    char *mode_ = NULL;
    char *fname_ = NULL;
}

phastaIO::~phastaIO () {
    if (mode_ != NULL) delete [] mode_;
    if (fname_ != NULL) delete [] fname_;
}

//
//  File Control
//

int phastaIO::openFile (const char *filename, const char *mode) {

    // Store the file name as a member variable:
    fileName_.append(filename);

    filePointer_=NULL ;
    fname_ = StringStripper( filename );
    mode_ = StringStripper( mode );

    if ( cscompare( mode_, "read" ) ) 
        filePointer_ = fopen( fname_, "rb" );
    else if( cscompare( mode_, "write" ) )
        filePointer_ = fopen( fname_ , "wb" );
    else if( cscompare( mode_, "append" ) )
        filePointer_ = fopen( fname_ , "ab" );

    //delete [] fname;
    //delete [] imode;

    if (filePointer_ == NULL) {
	perror("ERROR");
        fprintf(stderr,"ERROR opening file [%s].\n",fname_);
        return PHASTA_ERROR;
    }

    //fprintf(stdout,"fname_ : %s\n",fname_);
    //fprintf(stdout,"mode_ : %s\n",mode_);

    return PHASTA_OK;    

}

int phastaIO::closeFile() {
    if (cscompare(mode_, "write") || cscompare(mode_,"append")) {
      fflush(filePointer_);
    }
    fclose(filePointer_);
    return PHASTA_OK;
}

int phastaIO::rewindFile() {
     rewind(filePointer_);
     clearerr(filePointer_);
     return PHASTA_OK;
}


//
//  READ functions
//

int phastaIO::readHeader (const char* keyphrase,int* valueArray,
                          int  nItems,const char*  datatype,
                          const char*  iotype) {

   int i,skip_size,integer_value;
   int rewinded = 0;

   int lengthOfToken;
   originalPreColonTokenOnLine_.clear();
   originalPreColonTokenOnLine_.append("Not a token - reset at start of call to readHeader.");
   previousOriginalPreColonTokenOnLine_.erase();
   previousOriginalPreColonTokenOnLine_.append("Not a token - reset at start of call to readHeader.");

   isBinary( iotype );

   LastHeaderKey_[0] = '\0';

   while (rewinded < 2) {
     while (fgets(Line_, 1024, filePointer_) != NULL) {
         
         // ignore comment lines
         if (Line_[0] == '#') { 
             //fprintf(stdout,"COMMENT: %s",Line_);
             continue;
         }

         // ignore blank lines
         
         //char* chkblank = StringStripper(Line_);
         //fprintf(stdout,"chkblank: %i\n",strlen(chkblank));
         if (strlen(Line_) <= 1) {
             //fprintf(stdout,"Skipping blank line?\n");
             continue;
         }
         //fprintf(stdout,"Line_: %s",Line_);
         char* token = strtok ( Line_, ":" );
         
         // Get a copy of this key token, so that we can report it to terminal
         // in event of a read error...
         //
         // First, cycle back to save the previous token, too:
         previousOriginalPreColonTokenOnLine_ = originalPreColonTokenOnLine_;
         originalPreColonTokenOnLine_.erase();
         if (token != NULL)
         {
            originalPreColonTokenOnLine_.append(token);
         }
         else
         {
            originalPreColonTokenOnLine_.append("token was null");
         }

         //fprintf(stdout,"token: %s\n",token);
         if( cscompare( keyphrase , token ) ) {
            LastHeaderKey_[0] = '\0';
            sprintf(LastHeaderKey_,"%s",keyphrase); 
            token = strtok( NULL, " ,;<>" );
            skip_size = 0;
            skip_size = atoi( token );
            for( i=0;i < nItems && ( token = strtok( NULL," ,;<>") );i++) {
                valueArray[i] = atoi(token);
                //fprintf(stdout,"integer (%i): %i\n",i,atoi(token));
            }
            if ( i < nItems ) {
                fprintf(stderr,"Expected # of ints not recoverd from head\n");
                fprintf(stderr,"when looking for : %s\n", keyphrase);
                return PHASTA_ERROR;
            } else {
                return PHASTA_OK;
            }
         }
         // this really belongs when you open the file!
         if ( cscompare(token,"byteorder magic number") ) {
           if ( binary_format_ ) {
               fread(&integer_value, sizeof(int), 1, filePointer_);
               char junk;
               fread(&junk, sizeof(char), 1, filePointer_); /* reading the new line */
           } else {
               fgets(Line_,1024,filePointer_);
               sscanf(Line_,"%i",&integer_value);
           }
           if ( PHASTA_MAGIC_NUMBER != integer_value ) {
               //fprintf(stdout,"NOTE: Wrong_Endian found.\n");
             Wrong_Endian_ = true;
           }
           continue;
         }

         // skip to next header
         token = strtok( NULL, " ,;<>" );
         // Make sure that a token was actually found:
        if (token == NULL)
        {
            std::stringstream error;
            error << "EE: Failure whilst reading input file " << fileName_ << ", looking for key " << keyphrase << std::endl;
            std::cerr << "EE: Failure whilst reading input file. Failing looking for token: (probably not helpful!):" << std::endl;
            std::cerr << originalPreColonTokenOnLine_ << std::endl;
            std::cerr << "Perhaps more useful, the previous token that was found in the file was:" << std::endl;
            std::cerr << previousOriginalPreColonTokenOnLine_ << std::endl;
            throw std::runtime_error(error.str());
        }
         skip_size = atoi( token );
         if ( binary_format_ ) {
             fseek(filePointer_,skip_size,SEEK_CUR);
         } else {
            for( int gama=0; gama < skip_size; gama++ ) {
              fgets(Line_,1024,filePointer_);
            }
         }

     } // end inner while

     previousOriginalPreColonTokenOnLine_ = originalPreColonTokenOnLine_;
     originalPreColonTokenOnLine_.erase();
     originalPreColonTokenOnLine_.append("Not a token - the file was rewound.");

     rewind(filePointer_);
     if (ferror(filePointer_) != 0)
     {
        std::stringstream error;
        error << "EE: Failure after attempting to rewind " << fileName_ << ". ferror() detected error code " << ferror(filePointer_) << std::endl;
        throw std::runtime_error(error.str());
     }
     clearerr(filePointer_);
     rewinded++;

   }  // end outer while

   // std::cout << "Performing final (non-looping) rewind, having not found " << keyphrase << " in " << fileName_ << std::endl;
   // rewind(filePointer_);
   // if (ferror(filePointer_) != 0)
   // {
   //    std::stringstream error;
   //    error << "EE: Failure after attempting to rewind " << fileName_ << ". ferror() detected error code " << ferror(filePointer_) << std::endl;
   //    throw std::runtime_error(error.str());
   // }
   // clearerr(filePointer_);

   // std::cout << "final rewind allowed just occurred." << std::endl;

   // std::stringstream error;
   // error << "EE: Failure when attempting to read key " << keyphrase << " from file " << fileName_ << std::endl;
   // throw std::runtime_error(error.str());

   // This never gets returned if the error throw is still on the previous line!

   // WARNING - the above code for error throwing is something I added, then removed, when I discovered
   // that PHASTA_ERROR (returned here) sometimes means the simulation is OK to continue, and sometimes
   // means that it should terminate because the error is fatal. That's a real mess which should
   // be untangled some day, but it's a non-trivial job and there isn't time right now. CA 2015/1/9.
   return PHASTA_ERROR;

}


int phastaIO::readDataBlock ( const char* keyphrase,
                              void* valueArray,
                              int nItems,
                              const char*  datatype,
                              const char*  iotype ) {

    int n;
    int* valueArrayInt;	
    float* valueArrayFloat;
    double* valueArrayDouble;
    char junk;
   
    // check that the file has been opened
    if (filePointer_ == NULL) {
        fprintf(stderr,"No file associated with Descriptor \n");
        fprintf(stderr,"openfile_ function has to be called before \n");
        fprintf(stderr,"acessing the file\n");
        return PHASTA_ERROR;
    }
    
    if ( LastHeaderNotFound_ ) {
       fprintf(stderr,"ERROR:  last header was not found.\n");
       return PHASTA_ERROR;
    }

    // error check..
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.    

    if ( ! cscompare( LastHeaderKey_, keyphrase ) ) {
        fprintf(stderr,"ERROR: header not consistant with data block\n");
        fprintf(stderr,"Header: %s\n", LastHeaderKey_);
        fprintf(stderr,"DataBlock: %s\n", keyphrase);
        fprintf(stderr,"Please recheck read sequence\n");
        return PHASTA_ERROR;
    }

    size_t type_size = typeSize( datatype );
    isBinary( iotype );
    
    if ( binary_format_ ) {
        
        fread(valueArray,type_size,nItems,filePointer_ );
        fread(&junk, sizeof(char), 1, filePointer_); /* reading the new line */
        if ( Wrong_Endian_ ) SwapArrayByteOrder( valueArray, type_size, nItems );
        
    } else { 
        char junk;
        switch( type_of_data_ ) {
        case INT:
            valueArrayInt  = static_cast<int*>( valueArray );
            for( n=0; n < nItems ; n++ ) {
               int integer_value;
               Line_[0] = '\0';
               fgets(Line_,1024,filePointer_);
               if (sscanf(Line_,"%i",&integer_value) != 1) {
                   return PHASTA_ERROR;
               }
               valueArrayInt[n] = integer_value;
            }	
            break;
        case FLOAT:
            valueArrayFloat  = static_cast<float*>( valueArray );
            for( n=0; n < nItems ; n++ ) {
               float float_value;
               Line_[0] = '\0';
               fgets(Line_,1024,filePointer_);
               if (sscanf(Line_,"%f",&float_value) != 1) {
                   return PHASTA_ERROR;
               }
               valueArrayFloat[n] = float_value;
            }
            break;
        case DOUBLE:
            valueArrayDouble  = static_cast<double*>( valueArray );
            for( n=0; n < nItems ; n++ ) {
               double double_value;
               Line_[0] = '\0';
               fgets(Line_,1024,filePointer_);
               if (sscanf(Line_,"%lf",&double_value) != 1) {
                   return PHASTA_ERROR;
               }
               valueArrayDouble[n] = double_value;
            }
            break;
        }
    }
    
    return PHASTA_OK;
}

int phastaIO::readString(char* line) {
    while (fgets(line, 1024, filePointer_) != NULL) {
        return PHASTA_OK;
    }
    return PHASTA_ERROR;
}

//
//  WRITE functions
//

int phastaIO::writeString(const char* string) {
    fprintf(filePointer_,"%s",string);
    return PHASTA_OK;
}

int phastaIO::writeHeader (const char* keyphrase,
                           void* valueArray,
                           int nItems,
                           int ndataItems,
                           const char* datatype,
                           const char* iotype  ) {

    int* valueListInt;
   
    // check that the file has been opened
    if (filePointer_ == NULL) {
        fprintf(stderr,"No file associated with Descriptor \n");
        fprintf(stderr,"openfile_ function has to be called before \n");
        fprintf(stderr,"acessing the file\n");
        return PHASTA_ERROR;
    }

    LastHeaderKey_[0] = '\0';
    sprintf(LastHeaderKey_,"%s",keyphrase);

    DataSize_ = ndataItems;

    isBinary( iotype );
    size_t type_size = typeSize( datatype );
    
    int size_of_nextblock = 
        ( binary_format_ ) ? type_size*( ndataItems )+sizeof( char ) : ndataItems ;
    
    fprintf(filePointer_,"%s",keyphrase);
    fprintf(filePointer_," : < %i > ",size_of_nextblock);
    if( nItems > 0 ) {
        valueListInt = static_cast< int* > ( valueArray );
        for( int i = 0; i < nItems; i++ )
            fprintf(filePointer_,"%i ",valueListInt [i]);
    }
    fprintf(filePointer_,"\n");
    
    return PHASTA_OK ;
}


int phastaIO::writeDataBlock (const char* keyphrase,
                              void* valueArray,
                              int nItems,
                              const char* datatype,
                              const char* iotype ) {
    
    int n;

    // check that the file has been opened
    if (filePointer_ == NULL) {
        fprintf(stderr,"No file associated with Descriptor \n");
        fprintf(stderr,"openfile_ function has to be called before \n");
        fprintf(stderr,"acessing the file\n");
        return PHASTA_ERROR;
    }

    // error check..
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.    

    if ( ! cscompare( LastHeaderKey_, keyphrase ) ) {
        fprintf(stderr,"ERROR: header not consistant with data block\n");
        fprintf(stderr,"Header: %s\n", LastHeaderKey_);
        fprintf(stderr,"DataBlock: %s\n", keyphrase);
        fprintf(stderr,"Please recheck read sequence\n");
        return PHASTA_ERROR;
    }

    int* valueArrayInt;	
    float* valueArrayFloat;
    double* valueArrayDouble;
   
    int header_type = type_of_data_;
    size_t type_size=typeSize( datatype );

    if ( header_type != type_of_data_ ) {
        fprintf(stderr, "header and datablock differ on typeof data in the block for \n");
        fprintf(stderr, "keyphrase : %s\n", keyphrase);
        fprintf(stderr, "returning\n");
        return PHASTA_ERROR;
    }

    int nUnits = nItems;

    if ( nUnits != DataSize_ ) {
        fprintf(stderr, "header and datablock differ on the number of data items for\n");
        fprintf(stderr, "keyphrase : %s\n",keyphrase);
        fprintf(stderr, "returning\n");
        return PHASTA_ERROR;
    }
 
    isBinary( iotype );

    if ( binary_format_ ) {
        
      //size_t fwrite(const void* ptr, size_t size, size_t nobj, FILE* stream); 
        fwrite(static_cast< char* >( valueArray ),type_size, nUnits, filePointer_);
        fprintf(filePointer_,"\n");
        
    } else { 
        
        switch( type_of_data_ ) {
        case INT:
            
            valueArrayInt  = static_cast<int*>( valueArray );
            for( n=0; n < nUnits ; n++ ) 
                fprintf(filePointer_,"%i\n", valueArrayInt[n]);	
            break;

        case FLOAT:

            valueArrayFloat  = static_cast<float*>( valueArray );
            for( n=0; n < nUnits ; n++ ) 
                fprintf(filePointer_,"%f\n", valueArrayFloat[n]);	
            break;

        case DOUBLE:

            valueArrayDouble  = static_cast<double*>( valueArray );
            for( n=0; n < nUnits ; n++ ) 
                fprintf(filePointer_,"%lf\n", valueArrayDouble[n]);	
            break;
        }
    }	
    
    return PHASTA_OK;
}


//
//  Helper functions
//

char* phastaIO::StringStripper ( const char*  istring ) {
        char* fname;
        int namelength = strcspn( istring, " " );
        fname = new char [ namelength+1 ];
        strncpy( fname, istring , namelength );
        fname [ namelength ] = NULL;
        return fname;
}

int phastaIO::cscompare( const char* s1, const char* s2) {
        while( *s1 == ' ') s1++;
        while( *s2 == ' ') s2++;
        while( ( *s1 ) 
               && ( *s2 ) 
               && ( *s2 != '?')
               && ( tolower( *s1 )==tolower( *s2 ) ) ) {
            s1++;
            s2++;
            while( *s1 == ' ') s1++;
            while( *s2 == ' ') s2++;
        }
        if ( !( *s1 ) || ( *s1 == '?') ) return 1;
        else return 0;
}

void phastaIO::isBinary( const char* iotype ) {

        char* fname = StringStripper( iotype );
        if ( cscompare( fname, "binary" ) ) binary_format_ = true;
        else binary_format_ = false;
        delete [] fname;

}

size_t phastaIO::typeSize( const char* typestring ) {
        char* ts1 = StringStripper( typestring );
        if ( cscompare( ts1, "integer" ) ) {
            type_of_data_ = INT;
            delete [] ts1;
            return sizeof(int);
        } else if ( cscompare( ts1, "float" ) ) {
            type_of_data_ = FLOAT;
            delete [] ts1;
            return sizeof( float ); 
        } else if ( cscompare( ts1, "double" ) ) { 
            type_of_data_ = DOUBLE;
            delete [] ts1;
            return sizeof( double );
        } else { 
            delete [] ts1;
            fprintf(stderr,"unknown type\n");
            return 0;
        }
}

void phastaIO::SwapArrayByteOrder( void* array, int nbytes, int nItems ) {
        /* This swaps the byte order for the array of nItems each
           of size nbytes , This will be called only locally  */
        int i,j;
        unsigned char ucTmp;
        unsigned char* ucDst = (unsigned char*)array;
        
        for(i=0; i < nItems; i++) {
            for(j=0; j < (nbytes/2); j++)
                swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
            ucDst += nbytes;
        }
}

void openfile_( const char* filename, 
                const char* mode,
                int*  fileDescriptor ) {

    int i;

    // hard code allowable number of open files to 2048
    if (phastaIOfp == NULL) {
        phastaIOfp = new phastaIO* [2048];
        for (i = 0; i < 2048; i++) {
            phastaIOfp[i] = NULL;
        }
    }
    // skip 0 so it can be returned as an error code
    for (i = 1; i < 2048; i++) {
      if (phastaIOfp[i] == NULL) {
        phastaIOfp[i] = new phastaIO();
        if (phastaIOfp[i]->openFile(filename, mode) == PHASTA_ERROR) {
            delete phastaIOfp[i];
            phastaIOfp[i] = NULL;
            *fileDescriptor = 0;
            return;
        }
        *fileDescriptor = i;
        //fprintf(stdout,"file pointer (%i) opened\n",(*fileDescriptor));
        return;
      }
    }

    fprintf(stderr,"ERROR:  could not open file.\n");
    fprintf(stderr,"        maximum number of files exceeded.\n");
    exit(-1);
    return;

}

void closefile_( int* fileDescriptor, 
                 const char* mode ) {
    phastaIOfp[(*fileDescriptor)]->closeFile();
    delete phastaIOfp[(*fileDescriptor)];
    phastaIOfp[(*fileDescriptor)] = NULL;
    //fprintf(stdout,"file pointer (%i) closed\n",(*fileDescriptor));
    return;
}

int readheader_( int* fileDescriptor,
                  const char* keyphrase,
                  void* valueArray,
                  int*  nItems,
                  const char*  datatype,
                  const char*  iotype ) {
    int num = *nItems;
    return phastaIOfp[(*fileDescriptor)]->readHeader(keyphrase,(int*)valueArray,
                                            num,datatype,iotype);
}

void readdatablock_( int*  fileDescriptor,
                     const char* keyphrase,
                     void* valueArray,
                     int*  nItems,
                     const char*  datatype,
                     const char*  iotype ) {
    int num = *nItems;
    phastaIOfp[(*fileDescriptor)]->readDataBlock(keyphrase,
                       valueArray,num,datatype,iotype );
    return;
}

void writeheader_ (  int*  fileDescriptor,
                     const char* keyphrase,
                     void* valueArray,
                     int* nItems,
                     int* ndataItems,
                     const char* datatype,
                     const char* iotype  ) {

    int num = *nItems;
    int ndatanum = *ndataItems;
    phastaIOfp[(*fileDescriptor)]->writeHeader (keyphrase,
                     valueArray,num,ndatanum,datatype,iotype);
    return;
}

void writedatablock_( int* fileDescriptor,
                      const char* keyphrase,
                      void* valueArray,
                      int* nItems,
                      const char* datatype,
                      const char* iotype ) {
    int num = *nItems;
    phastaIOfp[(*fileDescriptor)]->writeDataBlock (keyphrase,
                             valueArray,num,datatype,iotype);
    return;
}

void writestring_( int* fileDescriptor,
                   const char* string ) {
    phastaIOfp[(*fileDescriptor)]->writeString(string);
}

void Gather_Headers( int* fileDescriptor, std::vector< std::string >& headers ) {

    char Line[1024];
    phastaIOfp[(*fileDescriptor)]->rewindFile();
     while ( phastaIOfp[(*fileDescriptor)]->readString(Line) == PHASTA_OK) {
        if ( Line[0] == '#' ) {
            headers.push_back( Line );
        } else { 
            break; 
        }
     }
    phastaIOfp[(*fileDescriptor)]->rewindFile();    
}

/*
int main (int argc, char **agrv) {

    int fd;
    openfile_("geombc.dat.0","read",&fd);
    int itwo=2;
    int* intfromfile = new int[3];
    readheader_(&fd,"co-ordinates?",intfromfile,&itwo,"double","binary?");
    int numnp=intfromfile[0];
    int nsd=intfromfile[1];
    double *xread = new double[numnp*nsd];
    int ixsiz=numnp*nsd;
    readdatablock_(&fd,"co-ordinates?",xread,&ixsiz,"double","binary?");
    return 0;
}
*/

void bzero_old(void* ptr, size_t sz) {
   char *cptr = (char*) ptr;
   for (int i=0; i < sz; i++) {
      cptr[i]=0;
   }
}
