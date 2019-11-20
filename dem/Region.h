//
// Created by maytee on 10/23/19.
//

#ifndef DEM_REGION_H
#define DEM_REGION_H

#include <vector>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <map>
#include <unordered_set>
#include "InputSettings.h"


// Data structure to hold a collection of grains. Can be used as a message between regions
// or as way to hold all neighbor grain data
struct GrainCollection {
  GrainCollection();
  void addGrainToCollection(double x_in, double y_in, double theta_in, double r_in, double vx_in, double vy_in, double omega_in);

  std::vector<double> q;
  std::vector<double> theta;
  std::vector<double> r;
  std::vector<double> v;
  std::vector<double> omega;
  std::vector<double> f;
  std::vector<int> unique_id;
  std::vector<std::vector<int>> neighbor_list;
  std::vector<std::vector<int>> contact_list;
  std::vector<std::vector<int>> contact_list_unique;
  std::vector<int> region_id;
  std::map<int,int> unique_to_local_map;
  int num_grains_in_collection;
};

struct ContactGrainGrain {
  int grain_idx_1;
  int grain_idx_2;
  std::vector<double> n{0.0,0.0};
  std::vector<double> contact_pt{0.0,0.0};
  double overlap;
  double s;
};

struct ContactGrainWall {
  int grain_idx_1;
  std::vector<double> n{0.0,0.0};
  std::vector<double> contact_pt{0.0,0.0};
  double overlap;
  double s;
};

class Region {
public:
  // Default Constructor
  Region();
  // Constructor to set up region bounds
  Region(std::vector<double> in_min, std::vector<double> in_max, int id, std::vector<bool> is_edge, double bin_size);
  // Constructor using input settings from input file
  Region(InputSettings& settings);


  // Simulation Functions ============================================================================
  // Generates Packing
  void generateRandomInitialPacking(double r_mean, double r_min, double r_max, int int_total_add);

  // Add new grain
  void addGrainToRegion(double x_in, double y_in, double theta_in, double r_in, double vx_in, double vy_in, double omega_in);
  // Add existing grain to region
  void addGrainToRegion(double x_in, double y_in, double theta_in, double r_in, double vx_in, double vy_in, double omega_in,
          double fx_in, double fy_in, std::vector<int> neighbor_list_in, std::vector<int> contact_list_in, int unique_id);

  // Add grain from region to collection
  void addGrainToCollection(GrainCollection& collection, int idx);

  // Add grain from one collection to another collection
  void addGrainFromCollectionToCollection(GrainCollection& collection_in, GrainCollection& collection_receive, int grain_idx);


  void removeGrain(int idx);

  // Create neighbor list
  void findNeighborsBruteForce(double delta);
  // Put grains into bins
  void rasterizeGrainsToBins(double delta);
  void findNeighborsFromBins(double delta);
  void findNeighborsFromBinsVect(double delta);

  // Find contacts
  void buildContactList();
  void buildContactListWithContactObjects();

  // Calculate and apply contact forces
  void calculateContactForces();
  void calculateContactForcesWithContactObjects();

  // Apply body forces
  void applyBodyForces();

  // Time integrate
  void forwardEuler(double delta_t);

  // Build map from unique ID to local ID
  void buildGrainIDMap();
  void buildGrainCollectionIDMap(GrainCollection& collection);

  // Communication Functions ============================================================================

  void clearCollectionData(GrainCollection& collection);

  // Builds a communication message of all grains that are within the cutoff distance away from the
  // boundary of whichever direction specified.
  // dir = 1 -> x, 2 -> y, 3 -> z.
  void buildCommunicationMessage(GrainCollection& message, int& dir, double& cutoff);

  void receiveCommunicationMessage(GrainCollection& in_message);

  void domainDecomposition(int num_regions, std::vector<Region>& region_list, std::vector<int>& domain_division);

//private:
  // Unique identifier for the region
  int m_region_id;

  // Bounds of region
  std::vector<double> m_region_min;
  std::vector<double> m_region_max;

  // Holds the list of grains that intersect each bin
  std::vector<std::vector<int>> m_bin_list;
  double m_bin_size;
  int m_num_bins_x;
  int m_num_bins_y;

  // Grain quantities for this region
  std::vector<double> m_q; // Position (qx1,qy1,qx2,qy2,..,qxn,qyn)
  std::vector<double> m_theta; // Rotational position (theta1,theta2,...,thetan)
  std::vector<double> m_r; // Radius (r1,r2,..,rn)
  std::vector<double> m_v; // Velocity (qx1,qy1,qx2,qy2,..,qxn,qyn)
  std::vector<double> m_omega; // Angular Velocity (omega1,omega2,...,omegan)
  std::vector<double> m_f; // Forces (qx1,qy1,qx2,qy2,..,qxn,qyn)
  std::vector<int> m_unique_id; // Unique ID of grain, for grain deletions and additions
  std::vector<std::vector<int>> m_neighbor_list; // List of neighbors for each grain
  std::vector<std::vector<int>> m_contact_list; // List of actual contacts for each grain
  std::vector<std::vector<int>> m_contact_list_unique; // List of actual contacts for each grain
  int m_num_grains; // Adjusts for grain addition and subtraction
  static int m_total_num_grains; // Lives across all regions
  static int m_id_tracker; // Lives across all regions so ID is unique across all regions
  std::vector<bool> m_is_edge; // Tells if it's an edge region or not. Order: [left,right,bottom,top]
  std::map<int,int> m_unique_to_local_map;

  // Quantities for surrounding grains
  GrainCollection m_surrounding_collection;


  std::map<std::pair<int,int>,ContactGrainGrain> m_contacts; // Key are the unique IDs
  std::map<std::pair<int,int>,ContactGrainGrain> m_contacts_cache; // Key are the unique IDs; from last timestep

  // Material parameters
  double k_n = 1000000.0;
  double k_t = 1000000.0;
  double eta_n = 75.55;
  double eta_t = 75.55;
  double m_mu = 0.5;
  double m_rho = 2600;
};


#endif //DEM_REGION_H
