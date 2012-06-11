/****************************************************************************** 

  (c) 2011-2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef _SCFIELD_H_
#define _SCFIELD_H_

namespace PHSOLVER {

class SCField {
   public:

   enum DomainType{ // http://msdn.microsoft.com/en-us/library/2dzy4k6e.aspx
      Volume,
      Surface,
      Scalar,
      DomainUndefined
   };
   enum ValType{
      Double,
      Integer,
      ValUndefined
   };

   SCField();

   SCField(char* name, DomainType domainType, int numUnits, int numVars, ValType valType);
   SCField(const char* pc_name, const char* pc_domainType, int numUnits, int numVars, const char* pc_valType);

   SCField(const SCField& inField);

   SCField& operator= (const SCField& inField);

   bool operator< (SCField inField) const;

   bool operator==(const SCField& inField) const;

   ~SCField();


   const char * GetName() const;
   DomainType GetDomainType() const;
   int GetNumVars() const;
   int GetNumUnits() const;
   ValType GetValType() const;

   int SetName(char* name);
   int SetDomainType(DomainType domainType);
   int SetNumVars(int numVars);
   int SetNumUnits(int numUnits);
   int SetValType(ValType valType);

   static const char* getDomainTypeName(DomainType dType);
   static const char* getValTypeName(ValType vType);

   private:
   char * pc_Name_; 
   DomainType e_DomainType_;  
   int i_NumVars_;
   int i_NumUnits_;
   ValType e_ValType_; 
   int ReadAndSetName_(const char* name);
   char* TerminateString_(const char* name);
};

}

#endif
