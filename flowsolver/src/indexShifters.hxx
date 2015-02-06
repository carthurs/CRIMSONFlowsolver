#ifndef INDEXSHIFTERS_HXX_
#define INDEXSHIFTERS_HXX_

inline int toZeroIndexing(const int oneIndexedValue)
{
	int zeroIndexedValue = oneIndexedValue - 1;
	return zeroIndexedValue;
}

inline int toOneIndexing(const int zeroIndexedValue)
{
	int oneIndexedValue = zeroIndexedValue + 1;
	return oneIndexedValue;
}

#endif