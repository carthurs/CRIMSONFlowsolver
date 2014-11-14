#include "dateTools.hxx"
#include <iostream>

void expiryDate::setExpiryDayOfMonth(int day_in)
{
	day = day_in;
}

void expiryDate::setExpiryMonthOfYear(int month_in)
{
	month = month_in;
}

void expiryDate::setExpiryYear(int year_in)
{
	year = year_in;
}

void expiryDate::checkWhetherExpiryDatePassed()
{
	writeExpiryDate();

	if (theDate.getYear() >= year &&
		theDate.getMonth() >= month &&
		theDate.getDay() > day)
		{
		    std::cout << "This copy of Simvascular Flowsolver has expired." << std::endl;
		    std::exit(1);
	    }
}

void expiryDate::writeExpiryDate()
{
	std::cout << "This copy of Simvascular Flowsolver expires on " << year << "/" << month << "/" << day << std::endl;
}

   // std::cout << "the year is " << theDate.now->tm_year <<std::endl;
   // std::cout << "the month is " << theDate.now->tm_mon <<std::endl;
   // std::cout << "the day is " << theDate.now->tm_mday <<std::endl;