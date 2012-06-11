/****************************************************************************** 

  (c) 2011-2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include "SCField.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>

namespace PHSOLVER {



SCField::SCField() {
  pc_Name_ = NULL;
  e_DomainType_ = SCField::DomainUndefined;
  i_NumVars_ = -1;
  i_NumUnits_ = -1;
  e_ValType_ = SCField::ValUndefined;
}

SCField::SCField(char* name, DomainType domainType, int numUnits,
                 int numVars, ValType valType) {
 
  ReadAndSetName_(name);
  e_DomainType_ = domainType;
  i_NumVars_ = numVars;
  i_NumUnits_ = numUnits;
  e_ValType_ = valType;
}

SCField::SCField(const char* pc_name, const char* pc_domainType, int numUnits,
                 int numVars, const char* pc_valType) {
  char* pc_cleanName = TerminateString_(pc_name);
  char* pc_cleanDomainType = TerminateString_(pc_domainType);
  char* pc_cleanValType = TerminateString_(pc_valType);
  
  ReadAndSetName_(pc_name);
  i_NumUnits_ = numUnits;
  i_NumVars_ = numVars;

  for ( int i_domainTypeIter = 0; i_domainTypeIter != SCField::DomainUndefined; i_domainTypeIter++ ) {
     if ( 0 == strcmp(pc_cleanDomainType,SCField::getDomainTypeName((SCField::DomainType) i_domainTypeIter)) ) {
       e_DomainType_ = (SCField::DomainType) i_domainTypeIter;
       break;
     }
  }
  SCField::ValType inputValType;
  for ( int i_valTypeIter = 0; i_valTypeIter != SCField::ValUndefined; i_valTypeIter++ ) {
     if ( 0 == strcmp(pc_cleanValType,SCField::getValTypeName((SCField::ValType) i_valTypeIter)) ) {
       e_ValType_ = (SCField::ValType) i_valTypeIter;
       break;
     }
  }

  free(pc_cleanName);
  free(pc_cleanDomainType);
  free(pc_cleanValType);
}


SCField::SCField(const SCField& inField) {
  ReadAndSetName_(inField.GetName());
  e_DomainType_ = inField.GetDomainType();
  i_NumVars_ = inField.GetNumVars();
  i_NumUnits_ = inField.GetNumUnits();
  e_ValType_ = inField.GetValType();
}

bool SCField::operator==(const SCField& inField) const {
  /*check that all vars are the same*/ 
  if ( 0 == strcmp(pc_Name_,inField.GetName())  &&
       e_DomainType_ == inField.GetDomainType() &&
       i_NumVars_ == inField.GetNumVars()       &&
       i_NumUnits_ == inField.GetNumUnits()       &&
       e_ValType_ == inField.GetValType()          ) {
    return true;
  }
  else {
    return false;
  }
}

SCField& SCField::operator= (const SCField& inField) {
  if ( this != &inField ) {
     ReadAndSetName_(inField.GetName());
     e_DomainType_ = inField.GetDomainType();
     i_NumVars_ = inField.GetNumVars();
     i_NumUnits_ = inField.GetNumUnits();
     e_ValType_ = inField.GetValType();
  }
  return *this;
}

bool SCField::operator< (SCField inField) const {
  if ( 0 > strcmp(pc_Name_, inField.GetName() ) ) {
    return true;
  }
  else {
    return false;
  }
}

const char* SCField::GetName() const { return pc_Name_; }

SCField::DomainType SCField::GetDomainType() const { return e_DomainType_; }

int SCField::GetNumVars() const { return i_NumVars_; }

int SCField::GetNumUnits() const { return i_NumUnits_; }

SCField::ValType SCField::GetValType() const { return e_ValType_; }

int SCField::SetName(char* name) {
  if (NULL != pc_Name_) {
    free(pc_Name_);
  }
  //ensure that the users input string is null terminated and has no whitespace at the end
  ReadAndSetName_(name);
  return 0;
}

int SCField::SetDomainType(DomainType domainType) { 
  e_DomainType_ = domainType; 
  return 0; 
}

int SCField::SetNumVars(int numVars) {
  i_NumVars_ = numVars; 
  return 0; 
}

int SCField::SetNumUnits(int numUnits) {
  i_NumUnits_ = numUnits; 
  return 0; 
}

int SCField::SetValType(ValType valType) {
  e_ValType_ = valType; 
  return 0; 
}

SCField::~SCField() {
  free(pc_Name_);
}

const char* SCField::getDomainTypeName(DomainType dType) {
   const char *pc_DomainTypeNames[] = {"Volume","Surface","Scalar"};
   return pc_DomainTypeNames[dType];
}
const char* SCField::getValTypeName(ValType vType) {
   const char* pc_ValTypeNames[] = {"Double","Integer"};
   return pc_ValTypeNames[vType];
}

char* SCField::TerminateString_(const char* pc_name) {
  #define MAX_STRING_LENGTH 2048
  char terminatedInputString[MAX_STRING_LENGTH];
  strncpy(terminatedInputString, pc_name, MAX_STRING_LENGTH-1);
  terminatedInputString[MAX_STRING_LENGTH-1] = 0;
  char *end = terminatedInputString + strlen(terminatedInputString) - 1;
  while (end >= terminatedInputString && 
         (isspace(*end) || iscntrl(*end)) ) {
    end--;
  }
  *(end + 1) = '\0';
  
  int inputStringLen = strlen(terminatedInputString);
  char * pc_terminatedString = (char*) calloc(inputStringLen+1, sizeof(char));
  strcpy(pc_terminatedString, terminatedInputString); 
  return pc_terminatedString;
}


int SCField::ReadAndSetName_(const char* name) {
  pc_Name_ = TerminateString_(name);
  return 0;
}

                 

} // END PHSOLVER NAMESPACE
