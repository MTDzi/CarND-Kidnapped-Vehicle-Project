/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"
#include <random>

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles;
	
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;

    // Random number generator engine
    std::default_random_engine gen;

    // Zero mean Gaussians
    std::normal_distribution<double> dist_x, dist_y, dist_theta;

    // Standard deviations for calculating Gaussian weights
    double std_lnd_mrk_x, var_x, std_lnd_mrk_y, var_y;
    double gaussian_normalizer;

    /**
     * one_particle_prediction Predicts the state for the next time step for a single particle
     * @param p The particle being modified
     * @param delta_t Time between time step t and t+1 in measurements [s]
     * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
     *   standard deviation of yaw [rad]]
     * @param velocity Velocity of car from t to t+1 [m/s]
     * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
     */
    void one_particle_prediction_(Particle & p, double dt, double v, double theta_dot);

    /**
	 * filter_landmarks_in_range_ Finds which map landmarks are in range
	 * @param sensor_range_2 Id est sensor_range * sensor_range
	 * @param x_p X-coordinate of a particle
     * @param y_p Y-coordinate of a particle
     * @param map Reference to the Map (which holds landmark observations)
	 */
    std::vector<LandmarkObs> filter_landmarks_in_range_(double sensor_range2, double x_p, double y_p, const Map & map);

    /**
	 * transform_observations_ Finds which map landmarks are in range
	 * @param x_p X-coordinate of a particle
     * @param y_p Y-coordinate of a particle
     * @param theta_p Angle at which the particle is headed
     * @param observations Reference to a vector of observations
	 */
    std::vector<LandmarkObs> transform_observations_(double x_p, double y_p, double theta_p, const std::vector<LandmarkObs> & observations);

    /**
	 * data_association_ Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
    void data_association_(std::vector<LandmarkObs> & in_map_obs, const std::vector<LandmarkObs> & landmark_obs);

public:
	
	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std_pos[], double std_landmark[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double velocity, double yaw_rate);

	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void update_weights(double sensor_range, const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and associations correctly connected
	 */
	Particle set_associations(Particle& particle, const std::vector<int>& associations,
		                     const std::vector<double>& sense_x, const std::vector<double>& sense_y);

	
	std::string get_associations(Particle best);
	std::string get_sense_X(Particle best);
	std::string get_sense_Y(Particle best);

	/**
	* initialized Returns whether particle filter is initialized yet or not.
	*/
	const bool initialized() const {
		return is_initialized;
	}
};



#endif /* PARTICLE_FILTER_H_ */
