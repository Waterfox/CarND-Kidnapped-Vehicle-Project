/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */


#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	normal_distribution<double> gaus_x_init(x, std[0]);
	normal_distribution<double> gaus_y_init(y, std[1]);
	normal_distribution<double> gaus_theta_init(theta, std[2]);
	default_random_engine gen2;

	//Number of Particles <IMPORTANT>
	num_particles = 10;
	particles.resize(num_particles);

	for (int i=0; i<num_particles; i++){
		particles[i].x = gaus_x_init(gen2);
		particles[i].y = gaus_y_init(gen2);
		// particles[i].y = y;
		// particles[i].x = x;
		// particles[i].theta = theta;
		particles[i].theta = gaus_theta_init(gen2);
		particles[i].weight = 1.0;
		weights.push_back(1.0f);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// We are passing sigma_pos (GPS noise) into this function as std_pos. Assume that this gets added in as
	// x,y noise and not as velocity / yaw rate noise per Lesson 7-20
	normal_distribution<double> x_noise(0.0, std_pos[0]);
	normal_distribution<double> y_noise(0.0, std_pos[1]);
	normal_distribution<double> theta_noise(0.0, std_pos[2]);
	default_random_engine gen3;


	for (int i=0; i<num_particles; i++){
		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double theta0 = particles[i].theta;

		if (yaw_rate > 0.000001) {
			particles[i].x = x0 + velocity/yaw_rate*(sin(theta0+yaw_rate*delta_t)-sin(theta0))+ x_noise(gen3);
			particles[i].y = y0 + velocity/yaw_rate*(cos(theta0)-cos(theta0+yaw_rate*delta_t))+ y_noise(gen3);
			particles[i].theta = theta0 + (yaw_rate*delta_t); // + theta_noise(gen3));
		}
		else {
			particles[i].x = x0 + velocity*cos(theta0)*delta_t + x_noise(gen3);
			particles[i].y = y0 + velocity*sin(theta0)*delta_t + y_noise(gen3);
			particles[i].theta = theta0 + (yaw_rate*delta_t) + theta_noise(gen3);
		}
	}
}

//SHOULD NOT WE BE COMPARING THESE TO THE MAP?
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	//iterate through each landmark in particle frame
	// TODO: Should check if there are zero predictions or observations
	for (int i=0; i<observations.size(); i++){
		double x_o = observations[i].x; //excess variables
		double y_o = observations[i].y;
		double min_dist = 999.9;
		double d_po = 999.9;
		//iterate through each observation in particle frame
		for (int j=0; j<predicted.size(); j++){

			double x_lm = predicted[j].x; //excess variables
			double y_lm = predicted[j].y;
			int id_lm = predicted[j].id;

			//if smallest distance, assign the lm ID to the obs
			d_po = dist(x_lm,y_lm,x_o,y_o);
			if (d_po < min_dist){
				min_dist = d_po;
				observations[i].id = id_lm;
			}
		}
		// cout << "observation " << i << " associated with LM " << observations[i].id << " dist=" << min_dist<<endl;
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
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	//iterate through each particle
	for (int i=0; i<num_particles; i++){

		//current particle's position
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;
		// cout <<"p_x: "<<p_x<<" p_y: "<<p_y<< " p_theta: " << p_theta << endl;


		//iterate through each observation and transform it to map frame
		// OBSERVATION CAR FRAME TO MAP FRAME
		std::vector<LandmarkObs> observations_map;

		for (int j=0; j<observations.size(); j++){

		  double x_o = observations[j].x;
		  double y_o = observations[j].y;
		  int id = observations[j].id; // shouldn't need this

		  double x_m;
		  double y_m;

		  //rotate by theta and translate
		  x_m = x_o*cos(p_theta) + y_o*-1.0*sin(p_theta) + p_x;
		  y_m = x_o*1.0*sin(p_theta) + y_o*cos(p_theta) + p_y; //NEGATIVE SIN FOR Y AXIS DOWN?

			LandmarkObs obs_m;
			obs_m.x = x_m;
			obs_m.y = y_m;
			obs_m.id = id;
			observations_map.push_back(obs_m);

		}

		// predicted map landmarks
		std::vector<LandmarkObs> lms_predicted;
		for (int k=0; k<map_landmarks.landmark_list.size(); k++) {

			if(dist(p_x, p_y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f) < sensor_range) {
				LandmarkObs lm_pred;
				lm_pred.id = map_landmarks.landmark_list[k].id_i;
				lm_pred.x = map_landmarks.landmark_list[k].x_f;
				lm_pred.y = map_landmarks.landmark_list[k].y_f;
				lms_predicted.push_back(lm_pred);
			}
		}

		dataAssociation(lms_predicted,observations_map);

		//calculate the weight of each particle
		double sigma_x = std_landmark[0];
		double sigma_y = std_landmark[1];
		double P = 1.0; // particles probability

		for (int k=0; k<observations_map.size(); k++){
			double x = observations_map[k].x;
			double y = observations_map[k].y;
			int id = observations_map[k].id; // associated LM ID
			double mu_x = 999.0;
			double mu_y = 999.0;
			// cout << "observation " << k << " associated with id " << id << endl;

			//search and match for the corresponding landmark ID
			for (int j=0; j<lms_predicted.size(); j++){
				// cout << "L id" << landmarks_p[j].id << endl;
				if (lms_predicted[j].id == id){
					mu_x = lms_predicted[j].x;
					mu_y = lms_predicted[j].y;
					// cout << "ASSOCIATED" << endl; //debug
					// cout << mu_x << " <x> " << x <<endl; //debug
					// cout << mu_y << " <y> " << y <<endl; //debug
					break;
				}
				//cout << "NO ASSOCIATION" << endl;
			}
			//update the total probability
			double p = 1.0/2.0/M_PI/sigma_x/sigma_y*exp(-1.0*(pow(x-mu_x,2)/(2.0*sigma_x)+ pow(y-mu_y,2)/(2.0*sigma_y)));
			// cout << "p= " << p << endl;
			P = P*p;
		}
		particles[i].weight = P;
		weights[i] = P;
		// cout << "P[" << i <<"]: " << P << endl;
	}

	//normalize weights 0..1 -- Don't need this
	// //find the max and min weights
	// double max_weight = 0.0;
	// double min_weight = 1.0;
	// for (int i=0; i<num_particles; i++){
	// 	if (particles[i].weight > max_weight){
	// 		max_weight = particles[i].weight;
	// 	}
	// 	if (particles[i].weight < min_weight){
	// 		min_weight = particles[i].weight;
	// 	}
	// }
	//
	// for (int i=0; i<num_particles; i++){
	// 	particles[i].weight = (particles[i].weight - min_weight)/(max_weight - min_weight);
	// 	// cout << "Pn[" << i <<"]: " << particles[i].weight << endl;
	// }


	// cout << "max_weight= " << max_weight << endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<> discdist(weights.begin(), weights.end());

	std::vector<Particle> resampled_particles;
	for (int i = 0; i < num_particles; ++i)
			resampled_particles.push_back(particles[discdist(gen)]);

	particles = resampled_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
