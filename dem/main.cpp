#include <iostream>
#include <vector>
#include "outputWriter.h"
#include "Region.h"
#include <random>
#include "chrono"
#include "omp.h"

int main(int argc, char *argv[]) {

  double r_mean = 0.05;
  double bin_size = r_mean*2.0*3.0;

  if(argc < 2) {
    std::cout<<"Using 3 grain d as bin size"<<std::endl;
  }
  else {
    bin_size = r_mean*2.0*std::atof(argv[1]);
    std::cout<<"Bin size of "<<std::atof(argv[1])<<" grain d"<<std::endl;
  }

  srand(100);
  std::cout << "Hello, World!" << std::endl;

  // Generate initial region
  std::vector<double> global_region_min{0.0,0.0};
  std::vector<double> global_region_max{24.0,16.0};
  int num_total_add = 32000;
  std::vector<bool> is_edge{1,1,1,1};
  Region init_region(global_region_min,global_region_max,0, is_edge,bin_size);
  init_region.generateRandomInitialPacking(r_mean,num_total_add);
  int num_total_init = init_region.m_num_grains;

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
  std::vector<int> domain_division;
  init_region.domainDecomposition(num_threads,region_list, domain_division);

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
  double t_f = 0.01;
  double dt = 0.00001;
  int num_regions = region_list.size();
  double delta = r_mean*1.21;
  int step_num = 0;
  int freq = 200;
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
    //std::cout<<"t: "<<t<<std::endl;
    step_num++;

    int right_dir = 1, left_dir = -1;
    double start_comms = omp_get_wtime();

    /*
    // Hardcode communication for now
    region_list[0].clearCollectionData(region_list[0].m_surrounding_collection);
    region_list[1].clearCollectionData(region_list[1].m_surrounding_collection);

     GrainCollection message1;
    region_list[1].buildCommunicationMessage(message1,left_dir,delta);
    region_list[0].receiveCommunicationMessage(message1);

    GrainCollection message2;
    region_list[0].buildCommunicationMessage(message2,right_dir,delta);
    region_list[1].receiveCommunicationMessage(message2);
    */

    int num_total_across_regions = 0;
    for(int ii = 0; ii < region_list.size(); ii++) {
     region_list[ii].clearCollectionData(region_list[ii].m_surrounding_collection);
     num_total_across_regions += region_list[ii].m_num_grains;
    }

    if(num_total_across_regions != num_total_init)
      std::cout << "LOST A GRAIN: "<<num_total_across_regions<<" vs "<<num_total_init<<std::endl;
    assert(num_total_across_regions == num_total_init);


    // Left-Right communication
    for(int x_dr = -1; x_dr < 2; x_dr += 2) {

      for (int y_div = 0; y_div < domain_division[1]; y_div++) {
        for (int x_div = 0; x_div < domain_division[0]; x_div++) {
          int region_num = x_div + y_div * domain_division[0];

          // If check left boundary and is left element, or check right boundary and check right element, skip
          if( (x_div == 0 && x_dr == -1) || x_div == domain_division[0]-1 && x_dr == 1) {
            continue;
          }
          //std::cout<<"LR Region "<<region_num<<" send message to "<<region_num+x_dr<<" "<<x_dr<<std::endl;
          GrainCollection message;
          region_list[region_num].buildCommunicationMessage(message,x_dr,delta);
          region_list[region_num+x_dr].receiveCommunicationMessage(message);
        }
      }
    }

    // Up-Down communication
    for(int y_dr = -1; y_dr < 2; y_dr += 2) {

      for (int y_div = 0; y_div < domain_division[1]; y_div++) {
        for (int x_div = 0; x_div < domain_division[0]; x_div++) {
          int region_num = x_div + y_div * domain_division[0];

          // If check bottom boundary and is bottom element, or check top boundary and check top element, skip
          if( (y_div == 0 && y_dr == -1) || y_div == domain_division[1]-1 && y_dr == 1) {
            continue;
          }
          std::cout<<"UD Region "<<region_num<<" send message to "<<region_num+y_dr*domain_division[0]<<" "<<y_dr<<std::endl;
          GrainCollection message;
          region_list[region_num].buildCommunicationMessage(message,y_dr,delta);
          region_list[region_num+y_dr*domain_division[0]].receiveCommunicationMessage(message);
        }
      }
    }
    double end_comms = omp_get_wtime();
    //std::cout<<"comm: "<<end_comms - start_comms<<std::endl;

    double start_map = omp_get_wtime();
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      region_list[region].buildGrainIDMap();
      region_list[region].buildGrainCollectionIDMap(region_list[region].m_surrounding_collection);
    }
    double end_map = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].rasterizeGrainsToBins(delta);
    }
    double end_raster = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].findNeighborsFromBins(delta);
    }
    double end_bins = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].buildContactList();
    }
    double end_contact = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].calculateContactForces();
    }
    double end_forces = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].applyBodyForces();
    }
    double end_body = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].forwardEuler(dt);
    }
    double end_forward = omp_get_wtime();

    std::cout<<end_map-start_comms<<","<<end_raster-end_map<<","<<end_bins-end_raster;
    std::cout<<","<<end_contact-end_bins<<","<<end_forces-end_contact;
    std::cout<<","<<end_body-end_forces<<","<<end_forward-end_body<<","<<end_forward-start_map<<std::endl;



    /*
    double start_for = omp_get_wtime();
    // Simulation in each region
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      //auto start_map = std::chrono::high_resolution_clock::now();
      region_list[region].buildGrainIDMap();
      region_list[region].buildGrainCollectionIDMap(region_list[region].m_surrounding_collection);
      //auto stop_map = std::chrono::high_resolution_clock::now();
      //auto duration_map = std::chrono::duration_cast<std::chrono::microseconds>(stop_map - start_map);
      //std::cout<<"map: "<<duration_map.count()<<std::endl;

      //auto start_rasterize = std::chrono::high_resolution_clock::now();
      region_list[region].rasterizeGrainsToBins(delta);
      //auto stop_rasterize = std::chrono::high_resolution_clock::now();
      //auto duration_raster = std::chrono::duration_cast<std::chrono::microseconds>(stop_rasterize - start_rasterize);
      //auto start_neighbor_bin = std::chrono::high_resolution_clock::now();
      //std::cout<<"grain raster: "<<duration_raster.count()<<std::endl;
      region_list[region].findNeighborsFromBins(delta);
      //auto stop_neighbor_bin = std::chrono::high_resolution_clock::now();
      //auto duration_neighbor_bin = std::chrono::duration_cast<std::chrono::microseconds>(stop_neighbor_bin - start_neighbor_bin);
      //std::cout<<"bin         : "<<duration_neighbor_bin.count()<<std::endl;
      region_list[region].buildContactList();
      region_list[region].calculateContactForces();
      region_list[region].applyBodyForces();
      region_list[region].forwardEuler(dt);
    }
    double end_for = omp_get_wtime();
    double num_contacts = 0;
    double num_in_bins = 0;

    //for(int region = 0; region < num_regions; region++){
    //  for(int ii = 0; ii < region_list[region].m_num_grains; ii++) {
    //    num_contacts += region_list[region].m_contact_list[ii].size();
    //  }
    //  for(int ii = 0; ii < region_list[region].m_num_bins_y*region_list[region].m_num_bins_x;ii++) {
    //    num_in_bins += region_list[region].m_bin_list[ii].size();
    //  }
    //}

    std::cout<<end_for - start_for<<",1";
    //std::cout<<","<<num_contacts<<","<<num_in_bins<<","<<end_comms-start_comms<<std::endl;
    std::cout<<","<<num_contacts<<","<<num_in_bins<<","<<end_comms-start_comms;
    std::cout<<","<<region_list[0].m_num_grains<<","<<region_list[0].m_surrounding_collection.num_grains_in_collection;
    std::cout<<","<<region_list[1].m_num_grains<<","<<region_list[1].m_surrounding_collection.num_grains_in_collection;
    std::cout<<std::endl;
    */


    bool verbose = false;

    if (verbose) {
      std::cout<<"Time: "<<t<<std::endl;
      for(int region = 0; region < num_regions; region++){
        for(int grain_num = 0; grain_num < region_list[region].m_num_grains; grain_num++){
          if(region_list[region].m_contact_list_unique[grain_num].size() > 0) {
            std::cout<<"Contact: "<<region_list[region].m_unique_id[grain_num]<<" with: ";
            for(int ii = 0; ii < region_list[region].m_contact_list[grain_num].size();ii++)
            {
              std::cout<<region_list[region].m_contact_list_unique[grain_num][ii]<<" ";
            }
            std::cout<<std::endl;

            std::cout<<"  Contact: "<<grain_num<<" with: ";
            for(int ii = 0; ii < region_list[region].m_contact_list[grain_num].size();ii++)
            {
              std::cout<<region_list[region].m_contact_list[grain_num][ii]<<" ";
            }
            std::cout<<std::endl;
          }
        }
      }
    }



    // Output to file
    if(step_num % skip == 0){
      //std::cout<<"output: "<<output_num<<" "<<t<<std::endl;
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