#include "CrimsonAugStatePart.h"

/*!
 *    This is the constructor for CrimsonAugStatePart.
 *    Currently, it sets the multiplying constant to 1.0
 *    and the name to an initial value.
 */
CrimsonAugStatePart::CrimsonAugStatePart()
:   premulconst_(1.0)
{
	name_ = "undefined";
}

/*!
 *    This is the destructor for CrimsonAugStatePart.
 *    Currently, nothing is explicitly done here.
 */
CrimsonAugStatePart::~CrimsonAugStatePart() {

}


/*!
      \param[in] const std::string& setname

      The only task that this function carries
      out currently is assigning the name of the
      augmented state part.
 */
void CrimsonAugStatePart::Initialize(const std::string& setname) {
	name_ = setname;
}


/*!
      \param[in] double* data_pointer

      This function adds a pointer to the
      internal vector of pointers.
 */
void CrimsonAugStatePart::addDataPointer(double* data_pointer) {
	pointers_.push_back(data_pointer);
}


/*!
      \param[in] bool val

      This function adds a boolean flag to
      the internal array of boolean flags denoting
      whether or not the associated quantity
      is included in the estimated portion
      of the augmented state.
 */
void CrimsonAugStatePart::addIsEstimated(bool val) {
	is_estimated_.push_back(val);
}


/*!
      \param[in] int ind

      This function assigns the index
      in CrimsonVerdandiModel::duplicated_state_
      that is used by CrimsonVerdandiModel::GetState()
      and CrimsonVerdandiModel::StateUpdated()
 */
void CrimsonAugStatePart::addDuplicatedStateIndex(int ind) {
	duplicated_state_indices_.push_back(ind);
}


/*!
      \param[in] int position
      \param[in] double* data_pointer

      This function assigns a pointer
      to a specific location in the internal array of pointers.
 */
void CrimsonAugStatePart::setDataPointer(int position, double* data_pointer) {
	pointers_[position] = data_pointer;
}


/*!
      \param[in] int position
      \param[in] double val

      This function assings a value to
      a location that is pointed to
      at the specified position
      in the internal array of pointers.
      Note that if the is_estimated_
      flag is 1, then, the actual
      assigned value is pow(2.0,val), i.e.
      the estimated quantity is re-parameterized.
 */
void CrimsonAugStatePart::setData(int position, double val) {
	is_estimated_[position] ? *(pointers_[position]) = pow(2.0,val)/premulconst_ : *(pointers_[position]) = val/premulconst_;
}



/*!
      \param[in] int position
      \param[in] bool val

      This function assigns the
      value of the is_estimated_ flag
      at the specified position.

 */
void CrimsonAugStatePart::setIsEstimated(int position, bool val) {
	is_estimated_[position] = val;
}


/*!
      \param[in] double val

      This function assigns the value
      for the constant that is multiplied
      to the output of getData().
 */
void CrimsonAugStatePart::setPremulConstant(double val) {
	premulconst_ = val;
}


/*!
      \param[in] int position

      This function returns the pointer
      at the specified position.
      Use this function to provide external access
      to the pointers_ array.
 */
double * CrimsonAugStatePart::getDataPointer(int position) {
	return pointers_[position];
}


/*!
      \param[in] int position

      This function returns the value pointed to
      at the specified position.
      If the associated is_estimated_ flag is true
      then the return value is actually log2 of the internal value.
 */
double CrimsonAugStatePart::getData(int position) {
	return is_estimated_[position] ? log2((*pointers_[position])*premulconst_) : (*pointers_[position])*premulconst_;
}


/*!
 *    This function returns
 *    the value of the multiplier constant
 */
double CrimsonAugStatePart::getPremulConstant() {
	return premulconst_;
}


/*!
      \param[in] int position

      This returns the value of the is_estimated_
      flag at a specified position
 */
bool CrimsonAugStatePart::getIsEstimated(int position) {
	return is_estimated_[position];
}


/*!
      \param[in] int position

      This returns the duplicated state index
      at a specified position
 */
int CrimsonAugStatePart::getDuplicatedStateIndex(int position) const {
	return duplicated_state_indices_[position];
}


/*!
 *    \return: std::size_t
 *
 *    This function returns the size of the
 *    internal pointers_ vector
 */
std::size_t CrimsonAugStatePart::getSize() const {
	return pointers_.size();
}


/*!
 *    \return: unsigned int
 *
 *    This function returns the total
 *    number of estimated variables in
 *    the augmented state component
 */
unsigned int CrimsonAugStatePart::getNumEstimated() {
	unsigned int sum_of_elems=0;
	for(std::vector<bool>::iterator j=is_estimated_.begin();j!=is_estimated_.end();++j)
	    sum_of_elems += *j;
	return sum_of_elems;
}


/*!
 *    \return: int
 *
 *    Returns the index of the first estimated
 *    variable
 *    in the augmented state component
 */
int CrimsonAugStatePart::getFirstEstimated() {
	std::vector<bool>::iterator j;
	int k;
	for(j=is_estimated_.begin(), k = 0;j!=is_estimated_.end();j++,k++) {
		if (*j) {
			return k;
		}
	}
	return -1;
}


/*!
 *    \return: std::string
 *
 *    Returns the name
 *    of the augmented state component
 */
std::string CrimsonAugStatePart::getName() const {
    return name_;
}


/*!
 *    This function resets the augmented state component
 *    by clearing the vectors and resetting
 *    the multiplying constant to 1.0
 *
 */
void CrimsonAugStatePart::Clear() {
	pointers_.clear();
	is_estimated_.clear();
	premulconst_ = 1.0;
}

