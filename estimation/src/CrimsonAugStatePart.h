#ifndef CrimsonAugStatePart_HXX
#define CrimsonAugStatePart_HXX

#include <string>
#include <vector>
#include <cmath>

//! This class organizes the pointers to the PHASTA arrays
//! into augmented state components
class CrimsonAugStatePart
{
public:

	//! Constructor.
	CrimsonAugStatePart();

	//! Destructor
	~CrimsonAugStatePart();

	//! Initialize with a name
	void Initialize(const std::string& setname);

	//! Adds a pointer to the array of pointers
	void addDataPointer(double* data_pointer);

	//! Adds to the array of boolean flags for the estimated variables
	void addIsEstimated(bool val);

	//! Sets the index in the duplicated state vector for a specific variable
    void addDuplicatedStateIndex(int ind);

	//! Sets a specific pointer in the array of pointers
	void setDataPointer(int position, double* data_pointer);

	//! Sets the data for a specific variable.
	void setData(int position, double val);

	//! Sets the estimated-variable-flag for a specific variable
	void setIsEstimated(int position, bool val);

	//! Sets the constant that the output of getData is multiplied with
	void setPremulConstant(double val);

	//! Returns the data pointer to a specific variable
	double * getDataPointer(int position);

	//! Returns the value of a specific variable
	double getData(int position);

	//! Returns the constant that the output of getData is multiplied with
	double getPremulConstant();

	//! Returns the estimated-variable-flag for a specific variable
	bool getIsEstimated(int position);

	//! Returns the index in the duplicated state vector for a specific variable
	int getDuplicatedStateIndex(int position) const;

	//! Returns the number of variables
	std::size_t getSize() const;

	//! Returns the number of estimated variables
	unsigned int getNumEstimated();

	//! Returns the index of the first estimated variable
    int getFirstEstimated();

    //! Returns the name
	std::string getName() const;

	//! Empty all the arrays and sets the multiplying constant to 1.0
	void Clear();

protected:

	//! Vector of pointers to PHASTA Fortran data
	std::vector<double*> pointers_;

	//! Vector of flags denoting if quantity is estimated
	std::vector<bool> is_estimated_;

	//! Vector of indices in the duplicated state
	std::vector<int> duplicated_state_indices_;

	//! Name of the augmented state component
	std::string name_;

	//! The constant that multiplies the
	double premulconst_;
};

#endif
