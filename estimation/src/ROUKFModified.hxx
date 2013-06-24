#ifndef ROUKFMODIFIED_HXX

#include "seldon/vector/VectorCollection.hxx"

//! This class inherits from the ROUKF class
template<class T, class Model, class ObservationManager>
class ROUKFModified: public Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager> {

public:

	//! Type of a row of the background error variance.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::model_state_error_variance_row
			model_state_error_variance_row;

	//! Type of the model state vector.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::model_state model_state;

	//! Type of the model/observation crossed matrix.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::matrix_state_observation
			matrix_state_observation;

	//! Type of the background error variance.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::model_state_error_variance
			model_state_error_variance;

	//! Type of the tangent linear model.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::model_tangent_linear_operator
			model_tangent_linear_operator;

	//! Type of the tangent linear observation operator.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>
			::observation_tangent_linear_operator observation_tangent_linear_operator;

	//! Type of the observation error variance.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>
			::observation_error_variance observation_error_variance;

	//! Type of a row of the tangent linear observation operator.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::observation_tangent_linear_operator_row
			observation_tangent_linear_operator_row;

	//! Type of the observation vector.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::observation
			observation;

	//! Type of the sigma point vector.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::sigma_point sigma_point;

	//! Type of the state vector collection.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::state_collection state_collection;

	//! Type of the observation vector collection.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::observation_collection observation_collection;

	//! Type of the sigma point matrix.
	typedef typename Verdandi::ReducedOrderUnscentedKalmanFilter<T,Model,ObservationManager>::sigma_point_matrix
	sigma_point_matrix;

	typedef typename Model::reduced_state_error_variance reduced_state_error_variance;

protected:

	int Nstate_local_;

public:

	ROUKFModified();
	~ROUKFModified();

	void GetReducedStateErrorVariance(reduced_state_error_variance &Preduced,
			int state_reduced_start_local);

	void SaveStateErrorVariance();

	void Analyze();

	void Forward();

};


#define ROUKFMODIFIED_HXX
#endif
