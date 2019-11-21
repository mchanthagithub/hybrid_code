//
// Created by maytee on 11/19/19.
//

#ifndef DEM_INPUTSETTINGS_H
#define DEM_INPUTSETTINGS_H

#include "vector"

struct InputSettings {
  std::vector<std::vector<double> > list_of_region_bounds;
  std::vector<std::vector<double> > list_radii_dist_params; // [mean, -%mean, +%mean, num]
  std::vector<std::vector<double> > list_of_individual_grains;
  std::vector<std::vector<double> > list_of_walls;

  double t_f = -1.0;
  double dt = -1.0;
  double k_n = -1.0;
  double k_t = -1.0;
  double eta_n = -1.0;
  double eta_t = -1.0;
  double m_mu = -1.0;
  double m_rho = -1.0;
  double bin_size = -1.0;
  int freq = -1;
  bool out_forces = 0;
};


#endif //DEM_INPUTSETTINGS_H
