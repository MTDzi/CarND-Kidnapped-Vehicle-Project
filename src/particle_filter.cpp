/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std_pos[], double std_landmark[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 10;

    dist_x = normal_distribution<double>(0, std_pos[0]);
    dist_y = normal_distribution<double>(0, std_pos[1]);
    dist_theta = normal_distribution<double>(0, std_pos[2]);
    std_lnd_mrk_x = std_landmark[0];
    std_lnd_mrk_y = std_landmark[1];
    var_x = std_lnd_mrk_x * std_lnd_mrk_x;
    var_y = std_lnd_mrk_y * std_lnd_mrk_y;
    gaussian_normalizer = 0.5 / M_PI / std_lnd_mrk_x / std_lnd_mrk_y;

    for (int i=0; i<num_particles; i++) {
        Particle one_particle = Particle();
        one_particle.id = i;
        one_particle.x = x + dist_x(gen);
        one_particle.y = y + dist_y(gen);
        one_particle.theta = theta + dist_theta(gen);
        one_particle.weight = 1.0;
        particles.push_back(one_particle);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double dt, double velocity, double theta_dot) {
	// TODO: Add measurements to each particle and add random Gaussian noise.

    for( auto & p: particles )
        one_particle_prediction_(p, dt, velocity, theta_dot);
}

void ParticleFilter::one_particle_prediction_(Particle & p, double dt, double velocity, double theta_dot) {

    if( std::fabs(theta_dot) > 0.0001 ) {
        double theta_inc = p.theta + dt * theta_dot;
        double v_by_theta_dot = velocity / theta_dot;

        p.x += v_by_theta_dot * (sin(theta_inc) - sin(p.theta));
        p.y += v_by_theta_dot * (cos(p.theta) - cos(theta_inc));
        p.theta = theta_inc;

    } else {
        double v_times_dt = velocity * dt;
        p.x += v_times_dt * cos(p.theta);
        p.y += v_times_dt * sin(p.theta);
//        theta stays the same, so there's no point in writing: p.theta = p.theta;
    }

    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);
}

void ParticleFilter::update_weights(double sensor_range, const std::vector<LandmarkObs> & observations, const Map & map) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node9 9.html

    double sensor_range2 = sensor_range * sensor_range;

    for( auto & part : particles ) {

        double x_p = part.x;
        double y_p = part.y;
        double theta_p = part.theta;

        // Leave only those Map landmarks that are within range (so it runs faster)
        std::vector<LandmarkObs> in_range_landmarks;
        in_range_landmarks = filter_landmarks_in_range_(sensor_range2, x_p, y_p, map);

        // Transform observations into the Map's coordinate system
        std::vector<LandmarkObs> transformed_observations;
        transformed_observations = transform_observations_(x_p, y_p, theta_p, observations);

        // Associate observations with landmarks on the map (set the id attribute for each transformed_observation)
        data_association_(transformed_observations, in_range_landmarks);

        // In this loop this particle's weight is being actually calculated
        part.weight = 1.0;

        for( auto & obs : transformed_observations ) {

            int associated_map_lnd_mrk = obs.id;
            double x1 = obs.x;
            double y1 = obs.y;

            double x2 = NAN;
            double y2 = NAN;

            // Find the associated map landmark, and extract it's X- and Y-coordinates
            for( auto & map_lnd_mrk : in_range_landmarks ) {

                if( map_lnd_mrk.id != associated_map_lnd_mrk )
                    continue;

                x2 = map_lnd_mrk.x;
                y2 = map_lnd_mrk.y;
                break;
            }
            if( std::isnan(x2) || std::isnan(y2) )
                std::cout << "Could not find an associated map landmark -- something went wrong";

            // Calculate the weight with a 2-D Gaussian
            double x_comp = pow(x1 - x2, 2) / 2 / var_x;
            double y_comp = pow(y1 - y2, 2) / 2 / var_y;
            double w = gaussian_normalizer * exp( -(x_comp + y_comp) );

            // And finally, modify this particle's weight
            part.weight *= w;
        }
    }
}

std::vector<LandmarkObs> ParticleFilter::filter_landmarks_in_range_(double sensor_range2, double x_p, double y_p, const Map & map) {

    std::vector<LandmarkObs> in_range_landmarks;

    for( auto & map_lnd_mrk : map.landmark_list ) {

        int id = map_lnd_mrk.id_i;
        float x_m = map_lnd_mrk.x_f;
        float y_m = map_lnd_mrk.y_f;

        double x_diff = x_p - x_m;
        double y_diff = y_p - y_m;

        double dist = x_diff * x_diff + y_diff * y_diff;

        if( dist < sensor_range2 )
            in_range_landmarks.push_back(LandmarkObs(id, x_m, y_m));
    }

    return in_range_landmarks;
}

std::vector<LandmarkObs> ParticleFilter::transform_observations_(double x_p, double y_p, double theta_p, const std::vector<LandmarkObs> & observations) {

    std::vector<LandmarkObs> transformed_observations;

    for( auto & obs : observations ) {

        int id = obs.id;
        double x_o = obs.x;
        double y_o = obs.y;

        double x_in_map = x_p + x_o * cos(theta_p) - y_o * sin(theta_p);
        double y_in_map = y_p + x_o * sin(theta_p) + y_o * cos(theta_p);

        transformed_observations.push_back(LandmarkObs(id, x_in_map, y_in_map));
    }

    return transformed_observations;
}

void ParticleFilter::data_association_(std::vector<LandmarkObs> & in_map_obs, const std::vector<LandmarkObs> & landmark_obs) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    for( auto & obs : in_map_obs ) {

        double min_dist = numeric_limits<double>::max();
        int map_id = -1;

        for( auto & lnd_mrk : landmark_obs ) {

            double x_diff = obs.x - lnd_mrk.x;
            double y_diff = obs.y - lnd_mrk.y;

            double dist = x_diff * x_diff + y_diff * y_diff;

            if( dist < min_dist ) {
                min_dist = dist;
                map_id = lnd_mrk.id;
            }
        }

        obs.id = map_id;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    weights.clear();
    for( auto & part : particles )
        weights.push_back(part.weight);

    std::discrete_distribution<> particle_distribution = std::discrete_distribution<>(weights.begin(), weights.end());
    std::vector<Particle> new_particles;
    for(int i=0; i<num_particles; i++)
        new_particles.push_back(particles[particle_distribution(gen)]);

    particles = new_particles;
}

string ParticleFilter::get_associations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " ") );
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::get_sense_X(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::get_sense_Y(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " ") );
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
