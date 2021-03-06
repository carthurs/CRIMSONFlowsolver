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

void expiryDate::enableExpiryDate()
{
	thisBuildHasAnExpiryDate = true;
}

void expiryDate::checkWhetherExpiryDatePassed()
{
	if (thisBuildHasAnExpiryDate)
	{
		writeExpiryDate();

		if ((theDate.getYear() > year) ||
			(theDate.getYear() == year && theDate.getMonth() > month) ||
			(theDate.getYear() == year && theDate.getMonth() == month && theDate.getDay() > day))
			{
				if (mMpiRank==0)
				{
			    	std::cout << "This copy of CRIMSON Flowsolver has expired." << std::endl;
				}
				MPI_Barrier(MPI_COMM_WORLD);
			    std::exit(1);
		    }
	}
	else
	{
		if (mMpiRank==0)
		{
			std::cout << "This CRIMSON Flowsolver binary has no built-in expiry date." << std::endl;
		}
	}
}

void expiryDate::writeExpiryDate()
{
	std::cout << "This copy of CRIMSON Flowsolver expires on " << year << "/" << month << "/" << day << std::endl;
	std::cout << "The date today is " << theDate.getYear() << "/" << theDate.getMonth() << "/" << theDate.getDay() << std::endl;
}

   // std::cout << "the year is " << theDate.now->tm_year <<std::endl;
   // std::cout << "the month is " << theDate.now->tm_mon <<std::endl;
   // std::cout << "the day is " << theDate.now->tm_mday <<std::endl;