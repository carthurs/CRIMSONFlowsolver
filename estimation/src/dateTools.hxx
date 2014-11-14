#ifndef DATETOOLS_HXX_
#define DATETOOLS_HXX_

#include <ctime>

class date
{
public:
	date()
	{
	}

	int getDay()
	{
		return day;
	}

	int getMonth()
	{
		return month;
	}
	
	int getYear()
	{
		return year;
	}

protected:
	time_t time;
	int year = -1;
	int month = -1;
	int day = -1;
};

class currentDate : public date
{
public:
	struct tm* now;
	currentDate()
	{
		// Time now:
		time = std::time(0);
		now = localtime(&time);
		day = now->tm_mday;
		month = now->tm_mon + 1;
		year = now->tm_year + 1900;
	}

private:
	time_t time;
};

class expiryDate : public date
{
public:
	expiryDate()
	{
		theDate = currentDate();
	}

	// These friendly functions take ints as they'd appear on a human calendar
	// (i.e. they hide any zero-indexing from you).
	void setExpiryDayOfMonth(int day_in);
	void setExpiryMonthOfYear(int month_in);
	void setExpiryYear(int year_in);

	void checkWhetherExpiryDatePassed();
	void writeExpiryDate();
private:
	currentDate theDate;

};

#endif