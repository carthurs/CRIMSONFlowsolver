#ifndef ROUKFMODIFIED_HXX

#include "seldon/vector/VectorCollection.hxx"

//! This class inherits from the ROUKF class
template<class T, class Model, class ObservationManager>
class ROUKFModified: public Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager> {

public:

	typedef typename Model::reduced_state_error_variance reduced_state_error_variance;

	ROUKFModified();
	~ROUKFModified();

	void GetReducedStateErrorVariance(reduced_state_error_variance &Preduced);

};


#define ROUKFMODIFIED_HXX
#endif
