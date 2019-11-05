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

  // Generate initial region
  std::vector<double> global_region_min{0.0,0.0};
  std::vector<double> global_region_max{4.0,10.5};
  int num_total_add = 200;
  double r_mean = 0.05;
  Region init_region(global_region_min,global_region_max,0);
  init_region.generateRandomInitialPacking(r_mean,num_total_add);

  // Get number of threads for domain decomposition
  int num_threads = 1;
  #pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }
  std::cout<<"Num threads: "<<num_threads<<std::endl;
  
  // Domain decomposition
  std::vector<std::vector<double> > region_mins;
  std::vector<std::vector<double> > region_maxes;
  std::vector<Region> region_list;
  init_region.domainDecomposition(num_threads,region_list);

  writeGrainsVTU("init_grains.vtu",init_region.m_q,init_region.m_v,
                   init_region.m_r,init_region.m_unique_id,init_region.m_f);

  for(int ii = 0; ii < region_list.size();ii++){
    std::cout<<"Region "<<region_list[ii].m_region_id<<" ngrains: "<<region_list[ii].m_num_grains<<std::endl;
  }

  /*
  std::vector<double> region1_min{0.0,0.0};
  std::vector<double> region1_max{2.0,10.5};
  std::vector<double> region2_min{2.0,0.0};
  std::vector<double> region2_max{4.0,10.5};

  Region region1(region1_min,region1_max,0);
  Region region2(region2_min,region2_max,1);

  int num_regions_y = 1;
  int num_regions_x = 2;

  region1.generateRandomInitialPacking(r_mean,num_total_add);
  region2.generateRandomInitialPacking(r_mean,num_total_add);
  //region1.addGrainToRegion(1.95,2.0,0.05,1.0,0.0);
  //region2.addGrainToRegion(2.07,2.0,0.05,0.0,0.0);


  writeRegionVTU("region1.vtu",region1.m_region_min,region1.m_region_max);
  writeRegionVTU("region2.vtu",region2.m_region_min,region2.m_region_max);

  std::vector<Region*> region_list;
  region_list.push_back(&region1);
  region_list.push_back(&region2);
  */

  double t = 0.0;
  double t_f = 0.5;
  double dt = 0.000005;
  int num_regions = region_list.size();
  double delta = r_mean*3.0;
  int step_num = 0;
  int freq = 400;
  //freq = static_cast<int>(1.0/dt);
  int skip = static_cast<int>(1.0/dt)/freq;

  int output_num = 0;

  // Write initial configuration
  for(int region_num = 0; region_num < num_regions; region_num++) {
    std::string prefix = "grain";
    prefix.append(std::to_string(region_num));
    prefix.append("_");
    std::string output_region_grain_str = generateOutputFileStr(prefix,freq,t_f,dt,output_num);
    writeGrainsVTU(output_region_grain_str,region_list[region_num].m_q,region_list[region_num].m_v,
                   region_list[region_num].m_r,region_list[region_num].m_unique_id,region_list[region_num].m_f);

    std::string region_prefix = "region";
    region_prefix.append(std::to_string(region_num));
    region_prefix.append(".vtu");
    writeRegionVTU(region_prefix,region_list[region_num].m_region_min,region_list[region_num].m_region_max);

  }

  output_num++;
  while(t < t_f){
    std::cout<<"t: "<<t<<std::endl;
    step_num++;

    // Hardcode communication for now
    region_list[0].clearCollectionData(region_list[0].m_surrounding_collection);
    region_list[1].clearCollectionData(region_list[1].m_surrounding_collection);
    int right_dir = 1, left_dir = -1;
    auto start_comms = std::chrono::high_resolution_clock::now();

    GrainCollection message1;
    region_list[1].buildCommunicationMessage(message1,left_dir,delta);
    region_list[0].receiveCommunicationMessage(message1);

    GrainCollection message2;
    region_list[0].buildCommunicationMessage(message2,right_dir,delta);
    region_list[1].receiveCommunicationMessage(message2);
    auto stop_comms = std::chrono::high_resolution_clock::now();
    auto duration_comms = std::chrono::duration_cast<std::chrono::microseconds>(stop_comms - start_comms);
    std::cout<<"comms:        "<<duration_comms.count()<<std::endl;

    /*
    int total_own = region_list[0].m_num_grains + region_list[1].m_num_grains;
    std::cout<<"Total is: "<<total_own<<std::endl;
    std::cout<<"Region 0 own: "<<region_list[0].m_num_grains<<std::endl;
    for(int ii = 0; ii < region_list[0].m_num_grains; ii++)
      std::cout<<region_list[0].m_unique_id[ii]<<std::endl;

    std::cout<<"Region 0 neighbors: "<<region_list[0].m_surrounding_collection.num_grains_in_collection<<std::endl;
    for(int ii = 0; ii < region_list[0].m_surrounding_collection.num_grains_in_collection; ii++)
      std::cout<<region_list[0].m_surrounding_collection.unique_id[ii]<<std::endl;

   std::cout<<"Region 1 own: "<<region_list[1].m_num_grains<<std::endl;
    for(int ii = 0; ii < region_list[1].m_num_grains; ii++)
      std::cout<<region_list[1].m_unique_id[ii]<<std::endl;

    std::cout<<"Region 1 neighbors: "<<region_list[1].m_surrounding_collection.num_grains_in_collection<<std::endl;
    for(int ii = 0; ii < region_list[1].m_surrounding_collection.num_grains_in_collection; ii++)
      std::cout<<region_list[1].m_surrounding_collection.unique_id[ii]<<std::endl;
    */

    // Simulation in each region
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      auto start_rasterize = std::chrono::high_resolution_clock::now();
      region_list[region].rasterizeGrainsToBins(delta);
      auto stop_rasterize = std::chrono::high_resolution_clock::now();
      auto duration_raster = std::chrono::duration_cast<std::chrono::microseconds>(stop_rasterize - start_rasterize);
      auto start_neighbor_bin = std::chrono::high_resolution_clock::now();
      std::cout<<"grain raster: "<<duration_raster.count()<<std::endl;
      region_list[region].findNeighborsFromBins(delta);
      auto stop_neighbor_bin = std::chrono::high_resolution_clock::now();
      auto duration_neighbor_bin = std::chrono::duration_cast<std::chrono::microseconds>(stop_neighbor_bin - start_neighbor_bin);
      std::cout<<"bin         : "<<duration_neighbor_bin.count()<<std::endl;
      region_list[region].buildContactList();
      region_list[region].calculateContactForces();
      region_list[region].applyBodyForces();
      region_list[region].forwardEuler(dt);
    }

    // Output to file
    if(step_num % skip == 0){
      std::cout<<"output: "<<output_num<<" "<<t<<std::endl;
      for(int region_num = 0; region_num < num_regions; region_num++) {
        std::string prefix = "grain";
        prefix.append(std::to_string(region_num));
        prefix.append("_");
        std::string output_region_grain_str = generateOutputFileStr(prefix,freq,t_f,dt,output_num);
        writeGrainsVTU(output_region_grain_str,region_list[region_num].m_q,region_list[region_num].m_v,
                       region_list[region_num].m_r,region_list[region_num].m_unique_id,region_list[region_num].m_f);

      }
      output_num++;
    }
    t += dt;
  }


  return 0;
}