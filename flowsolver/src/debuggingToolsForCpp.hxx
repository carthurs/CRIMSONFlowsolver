#ifndef DEBUGGINGTOOLSFORCPP_HPP_
#define DEBUGGINGTOOLSFORCPP_HPP_

#include <cstdio>
#include <iostream>
#include <cstdarg>

// Call MAGICAL_DEBUG() on a line to get that linenumber printed to console when the line is reached.
#define MAGICAL_DEBUG() magicalDebuggingUnicorn.writeToConsoleHere(__LINE__,__FILE__);

// To use these, call MAGICAL_whatever_PRINTER(4, typeWhatever1, typeWhatever2, typeWhatever3, typeWhatever4)
// DO NOT MIX TYPES IN ONE CALL - DOUBLES ONLY IN MAGICAL_DOUBLE_PRINTER, INTS ONLY IN MAGICAL_INT_PRINTER
// -- IF YOU MIX, THEN THE WHOLE CALL BEHAVES IN AN UNDEFINED WAY.
#define MAGICAL_DOUBLE_PRINTER(...) magicalDebuggingUnicorn.printSomeVariables<double>(__LINE__, __FILE__, __VA_ARGS__)
#define MAGICAL_INT_PRINTER(...) magicalDebuggingUnicorn.printSomeVariables<int>(__LINE__, __FILE__, __VA_ARGS__)

// HOW TO USE:
//
// * Include the header for this file in the file you want to debug
// * declare the debugger as:
//   magicalDebuggingUnicorn magicalDebuggingUnicorn;
//   in the place where you want to use it
// * write MAGICAL_DEBUG(); where you want the output to be
// * or write MAGICAL_INT_PRINTER(numberOfInts,myInt1,myInt2,myInt3, ... , myInt_numberOfInts )
// * or MAGICAL_DOUBLE_PRINTER(numberOfDoubles,myDouble1, ... , myDouble_numberOfDoubles )


class magicalDebuggingUnicorn
{
public:
	
	// File line is obtained from preprocessor macros
	void writeToConsoleHere(int fileLine, const char* fileName)
	{
		outputIndex++;
		printf("\nThis is debug print number %i, on line %i of file %s.\n",outputIndex,fileLine,fileName);
	}

	// This function is polymorphic, but must have only one type of output !
	template <typename Type> void printSomeVariables(int fileLine, const char* fileName, const int numberOfVariablesToPrint, ...)
	{
		// initialise the list of variables of unknown length
		va_list ellipsisVariables;

		// tell the macro that the ellipses comes after firstVarToPrint
		va_start(ellipsisVariables,numberOfVariablesToPrint);

		// Write some bonus info to the console:
		this->writeToConsoleHere(fileLine,fileName);

		std::cout << "===> DEBUG PRINTER: ";
		// -1 to compensate for the fact we have explicitly set the firstVarToPrint, which we did to get polymorphism
		for (int ii=0; ii<numberOfVariablesToPrint; ii++)
		{
			std::cout << va_arg(ellipsisVariables,Type) << " ";
		}
		std::cout << std::endl;
		va_end(ellipsisVariables);
	}

private:

	static int outputIndex;
};


#endif