/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // Set the number of particles.

           num_particles = 108;

           // Set standard deviations for x, y, theta.

           double std_x = std[0];
           double std_y = std[1];
           double std_theta = std[2];

           // This line creates a normal (Gaussian) distribution for x - GPS position

           normal_distribution<double> dist_x(x,std_x);

           // Create normal distribution for y and theta - GPS position

           normal_distribution<double> dist_y(y,std_y);

           normal_distribution<double> dist_theta(theta,std_theta);

           default_random_engine gen;

           //Initailaize all particles to first position (based on estimate of x,y, theta and their 
           // uncertainties from GPS

           for (int i = 0; i< num_particles; ++i){

           	       double sample_x,sample_y, sample_theta;

           	       // sample from the normal distribution to generate particles

           	       sample_x = dist_x(gen);
           	       sample_y = dist_y(gen);
           	       sample_theta = dist_theta(gen);

           	       Particle particle;

           	       particle.id = i;
           	       particle.x = sample_x;
           	       particle.y = sample_y;
           	       particle.theta = sample_theta;
           	       //initialize the weight to 1
           	       particle.weight = 1.0;

           	       particles.push_back(particle);
           	       weights.push_back(particle.weight);
           }

           is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	    double std_x = std_pos[0];
	    double std_y = std_pos[1];
	    double std_theta = std_pos[2];

	    default_random_engine gen;

	    // Create normal noise for x, y, theta around zero. Noise is propotional to timestep delta_t

	    normal_distribution<double> dist_x(0, std_x * delta_t);
	    normal_distribution<double> dist_y(0, std_y * delta_t);
	    normal_distribution<double> dist_theta(0, std_theta * delta_t);

	    for (int i = 0; i < num_particles ; ++i) {

	    	Particle &particle = particles[i];

	    	//calculate x, y position and add noise

	    	if (fabs(yaw_rate) < 0.00001){
	    		cout << "Zero Yaw rate ";
	    		particle.x += velocity * delta_t * cos(particle.theta) + dist_x(gen);
	    		particle.y += velocity * delta_t * sin(particle.theta) + dist_y(gen);
	    		particle.theta = particle.theta + dist_theta(gen);

	    	}
	    	else {
	    		particle.x += (velocity/yaw_rate) * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta)) + dist_x(gen);
	    		particle.y += (velocity/yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t)) + dist_y(gen);

	    		particle.theta += yaw_rate * delta_t + dist_theta(gen);

	    	}
	    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < observations.size(); i++){

		//Get current observation
		LandmarkObs obs = observations[i];

		//initialize minimum distance to maximum possible

		double min_dist = numeric_limits<double>::max();

		//initialize id of landmark from map placeholder to be associated with the observation

		int map_id = -1;

		for (unsigned int j = 0; j < predicted.size(); j++){
			//get current prediction
			LandmarkObs pred = predicted[j];

			// get distance between observed and predicted landmanrk
			double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);

			// find the predicted landmark nearer to observed landmark

			if (cur_dist < min_dist){
				min_dist = cur_dist;
				map_id = pred.id;
			}
		}

		observations[i].id = map_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i; i < num_particles; i++) {

		// Get the particle x, y coordinates
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		// create a vector to hold map landmark locations predicted to be within sensor range of the particle

		vector<LandmarkObs> predictions;

		for (unsigned int j = 0; j< map_landmarks.landmark_list.size();j++){

			// get id and x,y coordinates
			float lm_x = map_landmarks.landmark_list[j].x_f;
			float lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;

			//Let us consider landmarks with in sensor range

			if (fabs(lm_x - p_x) <= sensor_range && fabs(lm_y-p_y) <= sensor_range){

				predictions.push_back(LandmarkObs{lm_id,lm_x,lm_y});
			}
		}

		//create a copy of list of observations transformed from vehicle coordinates to map coordinates

		vector<LandmarkObs> transformed_obs;
		for (unsigned int j =0; j < observations.size(); j++){
			double transformed_x = cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y + p_x;
			double transformed_y = sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y + p_y;
			transformed_obs.push_back(LandmarkObs{observations[j].id,transformed_x,transformed_y});
		}
        // perform data association for the prediction and transformed obseravtions

		dataAssociation(predictions, transformed_obs);
	
		//reinitialize weight

		particles[i].weight = 1.0;

		for (unsigned int j =0; j < transformed_obs.size(); j++){

			// placeholder for observation and associated prediction coordinates

			double obs_x,obs_y, pred_x,pred_y;

			obs_x = transformed_obs[j].x;
			obs_y = transformed_obs[j].y;

			int associated_prediction = transformed_obs[j].id;

			// get the x,y co-ordinates of the prediction associated with current observation

			for (unsigned int k = 0; k < predictions.size(); k++){
				if (predictions[k].id == associated_prediction){
					pred_x = predictions[k].x;
					pred_y = predictions[k].y;

				}
			}

			// calculate weight for this observation with multivariate Gaussian

			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			double obs_weight = (1/(2*M_PI*std_x*std_y)) * exp ( - (pow(pred_x - obs_x,2)/(2*pow(std_x,2)) +
				(pow(pred_y - obs_y,2)/2*pow(std_y,2)))) ;

			// multiply all the weights to get particle weight
			particles[i].weight *=obs_weight;
		}	
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	vector<Particle> new_particles;

	// Get the current weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++){
		weights.push_back(particles[i].weight);
	}

	// generate random starting index for resampleing wheel
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	auto index = uniintdist(gen);

	//get max weight
	double max_weight = *max_element(weights.begin(), weights.end());

	//uniform random distribution
	uniform_real_distribution<double> unirealdist(0.0, max_weight);

	double beta = 0.0;

	//spin the wheel!

	for (int i = 0; i < num_particles; i++){
		beta += unirealdist(gen) * 2.0;
		while (beta > weights[index]){
			beta -= weights[index];
			index = (index + 1) % num_particles;

		}
		new_particles.push_back(particles[index]);

	}
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
