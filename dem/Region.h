//
// Created by maytee on 10/23/19.
//

#ifndef DEM_REGION_H
#define DEM_REGION_H

#include <vector>
#include <iostream>
#include <cmath>

// Data structure to hold a collection of grains. Can be used as a message between regions
// or as way to hold all neighbor grain data
struct GrainCollection {
  std::vector<double> q;
  std::vector<double> r;
  std::vector<double> v;
  std::vector<double> f;
  std::vector<int> unique_id;
  std::vector<std::vector<int>> neighbor_list;
  std::vector<std::vector<int>> contact_list;
  std::vector<int> region_id;
  int num_grains_in_collection;
};

class Region {
public:
  // Constructor to set up region bounds
  Region(std::vector<double> in_min, std::vector<double> in_max, int id);
  // Default Constructor
  Region();

  // Simulation Functions ============================================================================
  // Generates Packing
  void generateRandomInitialPacking(double r_mean, int int_total_add);

  // Add new grain
  void addGrainToRegion(double x_in, double y_in, double r_in, double vx_in, double vy_in);
  // Add existing grain to region
  void addGrainToRegion(double x_in, double y_in, double r_in, double vx_in, double vy_in, double fx_in, double fy_in,
                              std::vector<int> neighbor_list_in, std::vector<int> contact_list_in, int unique_id);

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

  // Find contacts
  void buildContactList();

  // Calculate and apply contact forces
  void calculateContactForces();

  // Apply body forces
  void applyBodyForces();

  // Time integrate
  void forwardEuler(double delta_t);

  // Communication Functions ============================================================================

  void clearCollectionData(GrainCollection& collection);

  // Builds a communication message of all grains that are within the cutoff distance away from the
  // boundary of whichever direction specified.
  // dir = 1 -> x, 2 -> y, 3 -> z.
  void buildCommunicationMessage(GrainCollection& message, int& dir, double& cutoff);

  void receiveCommunicationMessage(GrainCollection& in_message);

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
  std::vector<double> m_r; // Radius (r1,r2,..,rn)
  std::vector<double> m_v; // Velocity (qx1,qy1,qx2,qy2,..,qxn,qyn)
  std::vector<double> m_f; // Forces (qx1,qy1,qx2,qy2,..,qxn,qyn)
  std::vector<int> m_unique_id; // Unique ID of grain, for grain deletions and additions
  std::vector<std::vector<int>> m_neighbor_list; // List of neighbors for each grain
  std::vector<std::vector<int>> m_contact_list; // List of actual contacts for each grain
  int m_num_grains; // Adjusts for grain addition and subtraction
  static int m_total_num_grains; // Lives across all regions
  static int m_id_tracker; // Lives across all regions so ID is unique across all regions

  // Quantities for neighboring grains
  GrainCollection m_neighbor_collection;

  // Material parameters
  double k = 100000.0;
  double eta = 4.0;
  double m_rho = 2600;
};


#endif //DEM_REGION_H
