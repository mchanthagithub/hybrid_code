#include <iostream>
#include <vector>
#include "outputWriter.h"
#include "Region.h"
#include <random>
#include "chrono"
#include "omp.h"

int main() {
  srand(100);
  std::cout << "Hello, World!" << std::endl;
  std::vector<double> global_region_min{0.0,0.0};
  std::vector<double> global_region_max{4.0,10.5};

  std::vector<double> region1_min{0.0,0.0};
  std::vector<double> region1_max{2.0,10.5};
  std::vector<double> region2_min{2.0,0.0};
  std::vector<double> region2_max{4.0,10.5};

  Region region1(region1_min,region1_max,0);
  Region region2(region2_min,region2_max,1);

  int num_regions_y = 1;
  int num_regions_x = 2;

  int num_total_add = 400;
  double r_mean = 0.05;
  region1.generateRandomInitialPacking(r_mean,num_total_add);
  region2.generateRandomInitialPacking(r_mean,num_total_add);
  //region1.addGrainToRegion(1.95,2.0,0.05,1.0,0.0);
  //region2.addGrainToRegion(2.07,2.0,0.05,0.0,0.0);


  writeRegionVTU("region1.vtu",region1.m_region_min,region1.m_region_max);
  writeRegionVTU("region2.vtu",region2.m_region_min,region2.m_region_max);

  std::vector<Region*> region_list;
  region_list.push_back(&region1);
  region_list.push_back(&region2);

  double t = 0.0;
  double t_f = 0.5;
  double dt = 0.00001;
  int num_regions = 2;
  double delta = r_mean*3.0;
  int step_num = 0;
  int freq = 400;
  int skip = static_cast<int>(1.0/dt)/freq;

  int output_num = 0;

  // Write initial configuration
  std::string output_region1_grain_str = generateOutputFileStr("grain1_",freq,t_f,dt,output_num);
  std::string output_region2_grain_str = generateOutputFileStr("grain2_",freq,t_f,dt,output_num);

  writeGrainsVTU(output_region1_grain_str,region1.m_q,region1.m_v,region1.m_r,region1.m_unique_id,region1.m_f);
  writeGrainsVTU(output_region2_grain_str,region2.m_q,region2.m_v,region2.m_r,region2.m_unique_id,region1.m_f);
  output_num++;
  bool brute = false;
  while(t < t_f){
    std::cout<<"t: "<<t<<std::endl;
    step_num++;


    // Hardcode communication for now
    region_list[0] -> clearCollectionData(region_list[0]->m_surrounding_collection);
    region_list[1] -> clearCollectionData(region_list[1]->m_surrounding_collection);
    int right_dir = 1, left_dir = -1;
    auto start_comms = std::chrono::high_resolution_clock::now();

    GrainCollection message1;
    region_list[1] -> buildCommunicationMessage(message1,left_dir,delta);
    region_list[0] -> receiveCommunicationMessage(message1);

    GrainCollection message2;
    region_list[0] -> buildCommunicationMessage(message2,right_dir,delta);
    region_list[1] -> receiveCommunicationMessage(message2);
    auto stop_comms = std::chrono::high_resolution_clock::now();
    auto duration_comms = std::chrono::duration_cast<std::chrono::microseconds>(stop_comms - start_comms);
    std::cout<<"comms:        "<<duration_comms.count()<<std::endl;


    // Simulation in each region
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      if(brute) {
        auto start_brute = std::chrono::high_resolution_clock::now();
        region_list[region]->findNeighborsBruteForce(delta);
        auto stop_brute = std::chrono::high_resolution_clock::now();
        auto duration_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_brute - start_brute);
        std::cout<<"brute: "<<duration_brute.count()<<std::endl;
      }
      else {
        auto start_rasterize = std::chrono::high_resolution_clock::now();
        region_list[region]->rasterizeGrainsToBins(delta);
        auto stop_rasterize = std::chrono::high_resolution_clock::now();
        auto duration_raster = std::chrono::duration_cast<std::chrono::microseconds>(stop_rasterize - start_rasterize);
        auto start_neighbor_bin = std::chrono::high_resolution_clock::now();
        std::cout<<"grain raster: "<<duration_raster.count()<<std::endl;
        region_list[region]->findNeighborsFromBins(delta);
        auto stop_neighbor_bin = std::chrono::high_resolution_clock::now();
        auto duration_neighbor_bin = std::chrono::duration_cast<std::chrono::microseconds>(stop_neighbor_bin - start_neighbor_bin);
        std::cout<<"bin         : "<<duration_neighbor_bin.count()<<std::endl;
      }
      region_list[region]->buildContactList();
      region_list[region]->calculateContactForces();
      region_list[region]->applyBodyForces();
      region_list[region]->forwardEuler(dt);
    }

    // Output to file
    if(step_num % skip == 0){
      std::cout<<"output: "<<output_num<<" "<<t<<std::endl;
      std::string output_region1_grain_str = generateOutputFileStr("grain1_",freq,t_f,dt,output_num);
      std::string output_region2_grain_str = generateOutputFileStr("grain2_",freq,t_f,dt,output_num);

      writeGrainsVTU(output_region1_grain_str,region1.m_q,region1.m_v,region1.m_r,region1.m_unique_id,region1.m_f);
      writeGrainsVTU(output_region2_grain_str,region2.m_q,region2.m_v,region2.m_r,region2.m_unique_id,region2.m_f);
      output_num++;
    }
    t += dt;
  }


  return 0;
}