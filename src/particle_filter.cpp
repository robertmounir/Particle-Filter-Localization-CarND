/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100; // TODO: Set the number of particles

  default_random_engine gen;
  double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

  // Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; ++i)
  {

    struct Particle p;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
  }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  default_random_engine gen;
  double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
                                  // Set standard deviations for x, y, and theta
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];
  double x;
  double y;
  double theta;

  normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);

  for (int i = 0; i < num_particles; ++i)
  {

     if (fabs(yaw_rate) < 0.00001) {  
      x=particles[i].x + velocity * delta_t * cos(particles[i].theta);
      y=particles[i].y + velocity * delta_t * sin(particles[i].theta);
      theta= particles[i].theta;
    }else
    {
    x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
    y = particles[i].y + (velocity / yaw_rate) * (-1 * cos(particles[i].theta + (yaw_rate * delta_t)) + cos(particles[i].theta));
    theta = particles[i].theta + delta_t * yaw_rate;
    }

    particles[i].x = dist_x(gen)+x;
    particles[i].y = dist_y(gen)+y;
    particles[i].theta = dist_theta(gen)+theta;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
/*
  vector<LandmarkObs> observations_p;
  for (int i = 0; i < predicted.size(); i++)
  {
    int id = predicted[i].id;
    int indx = 0;
    double dis = -1;
    for (int j = 0; j < observations.size(); j++)
    {
      if (observations[j].id = id)
      {
        double d = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
        if (dis == -1)
        {
          dis = d;
          indx = j;
        }
        else
        {
          if (d < dis)
          {
            dis = d;
            indx = j;
          }
        }
      }
    }
    observations_p.push_back(observations[indx]);
  }
  observations = observations_p;
  */
   for (unsigned int i = 0; i < observations.size(); i++) {
    
    // grab current observation
    LandmarkObs o = observations[i];

    // init minimum distance to maximum possible
    double min_dist = numeric_limits<double>::max();

    // init id of landmark from map placeholder to be associated with the observation
    int map_id = -1;
    
    for (unsigned int j = 0; j < predicted.size(); j++) {
      // grab current prediction
      LandmarkObs p = predicted[j];
      
      // get distance between current/predicted landmarks
      double cur_dist = dist(o.x, o.y, p.x, p.y);

      // find the predicted landmark nearest the current observed landmark
      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        map_id = p.id;
      }
    }

    // set the observation's id to the nearest predicted landmark's id
    observations[i].id = map_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (int i = 0; i < num_particles; ++i)
  {

    Particle p = particles[i];
    vector<LandmarkObs> predicted;
    vector<LandmarkObs> observations_p = observations;

    for (int x = 0; x < observations_p.size(); x++)
    {
      observations_p[x].x = p.x + cos(p.theta) * observations[x].x - sin(p.theta ) * observations[x].y;
      observations_p[x].y = p.y + sin(p.theta ) * observations[x].x + cos(p.theta ) * observations[x].y;
    }

    for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      if (dist(p.x, p.y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range)
      {
        LandmarkObs obj;
        obj.id = map_landmarks.landmark_list[j].id_i;
        obj.x = map_landmarks.landmark_list[j].x_f;
        obj.y = map_landmarks.landmark_list[j].y_f;
        predicted.push_back(obj);
      }
    }

    dataAssociation(predicted, observations_p);
    double sx = std_landmark[0];
    double sy = std_landmark[1];
     p.weight = 1;
    for (int j = 0; j < observations_p.size(); j++)
    {
     
       double o_x, o_y, pr_x, pr_y;
      o_x = observations_p[j].x;
      o_y = observations_p[j].y;

      int associated_prediction = observations_p[j].id;

      for (unsigned int k = 0; k < predicted.size(); k++) {
        if (predicted[k].id == associated_prediction) {
          pr_x = predicted[k].x;
          pr_y = predicted[k].y;
        }
      }




      double dx = (pr_x -o_x);
      double dy = (pr_y -o_y);
      p.weight *= (1 / (2 * M_PI * sx * sy)) * exp(-((dx * dx / 2 * sx * sx) + (dy * dy / 2 * sy * sy)));
    }

    particles[i].weight = p.weight;
  
  }



}

void ParticleFilter::resample()
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

 vector<Particle> new_particles;

  default_random_engine gen;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}