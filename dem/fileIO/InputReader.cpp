//
// Created by maytee on 11/14/19.
//

#include <vector>
#include "InputReader.h"

template<typename T> T getValueFromEqualityString(const char* value_char, std::string line) {
  std::string value = value_char;
  std::string delim = "=";
  if(line.find(value) != std::string::npos) {
    // Find everything after the equal sign
    std::string token = line.substr(line.find(delim)+1,std::string::npos);
    if(std::is_same<T,int>::value) {
      //std::cout<<value_char<<":"<<token<<std::endl;
      return std::stoi(token);
    } else if(std::is_same<T,double>::value) {
      //std::cout<<value_char<<" "<<token<<std::endl;
      return std::stod(token);
    }
  } else {
    std::cout<<"Unable to find setting "<<value<<" from input file.\nRead-in line is: "<<line<<std::endl;
  }
}

template<typename T> std::vector<T> getValuesFromList(std::string line) {
  std::string delim = ",";
  std::vector<T> data;
  size_t pos = 0;
  //std::cout<<"New data: ";
  int ctr = 0;
  while( (pos = line.find(delim)) != std::string::npos) {
    std::string token = line.substr(0, pos);
    //std::cout << token << " ";
    if (std::is_same<T, int>::value) {
      data.push_back(std::stoi(token));
    } else if (std::is_same<T, double>::value) {
      data.push_back(std::stod(token));
    }
    line.erase(0, pos + delim.length());
    ctr++;
  }
  // Get last value
  std::string token = line.substr(0, pos);
  //std::cout << token << " ";
  if (std::is_same<T, int>::value) {
    data.push_back(std::stoi(token));
  } else if (std::is_same<T, double>::value) {
    data.push_back(std::stod(token));
  }

  //std::cout<<std::endl;

  if(ctr == 0) {
    std::cout<<"Did not read in any list data. Line read is: "<<line<<std::endl;
    exit(0);
  }
  return data;
}


template<typename T> std::vector<T> extractData(std::ifstream& input_file, std::string field) {

}

template<typename T> void valueCheck(T value, std::string value_str) {
  if(value > 0.0) {
    std::cout<<" "<<value_str<<": "<<value<<std::endl;
  }
  else {
    std::cout<<" "<<value_str<<" not read in properly. Exiting"<<std::endl;
    exit(0);
  }
}

InputSettings readInputFile(std::string file_name) {
  InputSettings read_in_settings;
  std::string line;
  std::ifstream input_file (file_name);
  if(input_file.is_open()) {
    while(getline(input_file,line)) {
      if(line.find("#") != std::string::npos) {
      // Found a comment so skip
        continue;
      } else if(line.find("generate_rectangle=") != std::string::npos) {
        int num_generate = getValueFromEqualityString<int>("generate_rectangle",line);
        std::vector<std::vector<double> > region_bounds;
        for(int ii = 0; ii < num_generate; ii++) {
          getline(input_file,line);
          // Skip any comments
          while(line.find("#") != std::string::npos) {
            getline(input_file,line);
          }
          region_bounds.push_back(getValuesFromList<double>(line));
        }

        // Get radii distribution parameters
        std::vector<std::vector<double> > radii_distribution_params;
        for(int ii = 0; ii < num_generate; ii++) {
          getline(input_file,line);
          // Skip any comments
          while(line.find("#") != std::string::npos) {
            getline(input_file,line);
          }
          radii_distribution_params.push_back(getValuesFromList<double>(line));
        }

        if(region_bounds.size() != num_generate ) {
          std::cout<<"Number of rectangles does not match input number: "<<region_bounds.size()<<" vs "<<num_generate<<std::endl;
          exit(0);
        } else if (region_bounds.size() != radii_distribution_params.size()) {
          std::cout<<"Number of rectangles does not match number of distribution params: "<<region_bounds.size()<<" vs "<<radii_distribution_params.size()<<std::endl;
          exit(0);
        }
        read_in_settings.list_of_region_bounds = region_bounds;
        read_in_settings.list_radii_dist_params = radii_distribution_params;

      } else if(line.find("grain_list=") != std::string::npos) {
        int num_generate = getValueFromEqualityString<int>("grain_list",line);
        std::vector<double> grain_properties;
        for(int ii = 0; ii < num_generate; ii++) {
          getline(input_file,line);
          grain_properties = getValuesFromList<double>(line);
          if(grain_properties.size() != 7) {
            std::cout<<"Not enough grain properties given for a grain. Must provide: x,y,theta,r,vx,vy,omega"<<std::endl;
            exit(0);
          }
          read_in_settings.list_of_individual_grains.push_back(grain_properties);
        }
      } else if(line.find("walls=") != std::string::npos) {
        int num_walls = getValueFromEqualityString<int>("walls",line);
        std::vector<std::vector<double> > wall_list;
        for(int ii = 0; ii < num_walls; ii++) {
          getline(input_file,line);
          // Skip any comments
          while(line.find("#") != std::string::npos) {
            getline(input_file,line);
          }
          if(line.find("=") != std::string::npos) {
            std::cout<<"Found another property while going through wall list: "<<line<<std::endl;
            exit(0);
          }
          wall_list.push_back(getValuesFromList<double>(line));
        }

        if(wall_list.size() != num_walls ) {
          std::cout<<"Number of walls does not match input number: "<<wall_list.size()<<" vs "<<num_walls<<std::endl;
          exit(0);
        }
        read_in_settings.list_of_walls = wall_list;
      }else if (line.find("k_n=") != std::string::npos) {
        read_in_settings.k_n = getValueFromEqualityString<double>("k_n",line);
      }else if (line.find("k_t=") != std::string::npos) {
        read_in_settings.k_t = getValueFromEqualityString<double>("k_t",line);
      }else if (line.find("eta_n=") != std::string::npos) {
        read_in_settings.eta_n = getValueFromEqualityString<double>("eta_n",line);
      }else if (line.find("eta_t=") != std::string::npos) {
        read_in_settings.eta_t = getValueFromEqualityString<double>("eta_t",line);
      }else if (line.find("mu=") != std::string::npos) {
        read_in_settings.m_mu = getValueFromEqualityString<double>("mu",line);
      }else if (line.find("rho=") != std::string::npos) {
        read_in_settings.m_rho = getValueFromEqualityString<double>("rho",line);
      }else if (line.find("t_f=") != std::string::npos) {
        read_in_settings.t_f = getValueFromEqualityString<double>("t_f",line);
      }else if (line.find("dt=") != std::string::npos) {
        read_in_settings.dt = getValueFromEqualityString<double>("dt",line);
      }else if (line.find("bin_size=") != std::string::npos) {
        read_in_settings.bin_size = getValueFromEqualityString<double>("bin_size",line);
      }else if (line.find("out_freq=") != std::string::npos) {
        read_in_settings.freq = getValueFromEqualityString<int>("out_freq",line);
      }else if (line.find("out_forces=") != std::string::npos) {
        read_in_settings.out_forces = static_cast<bool>(getValueFromEqualityString<int>("out_forces",line));
      }
    }

  } else {
    std::cout<<"Could not open DEM input file "<<file_name<<std::endl;
    exit(0);
  }
  input_file.close();

  std::cout<<"Settings read in as: "<<std::endl;
  std::cout<<"Regions: "<<read_in_settings.list_of_region_bounds.size()<<std::endl;
  for(int ii = 0; ii < read_in_settings.list_of_region_bounds.size();ii++) {
    std::cout<<" "<<ii<<": ";
    for(int jj = 0; jj < read_in_settings.list_of_region_bounds[ii].size(); jj++) {
      std::cout<<read_in_settings.list_of_region_bounds[ii][jj];
      if(jj < read_in_settings.list_of_region_bounds[ii].size()-1)
        std::cout<<",";
    }
    std::cout<<std::endl;
  }
  std::cout<<"Number of individual grains listed: "<<read_in_settings.list_of_individual_grains.size()<<std::endl;

  std::cout<<"Walls: "<<read_in_settings.list_of_walls.size()<<std::endl;
  for(int ii = 0; ii < read_in_settings.list_of_walls.size();ii++) {
    std::cout<<" "<<ii<<": ";
    for(int jj = 0; jj < read_in_settings.list_of_walls[ii].size(); jj++) {
      std::cout<<read_in_settings.list_of_walls[ii][jj];
      if(jj < read_in_settings.list_of_walls[ii].size()-1)
        std::cout<<",";
    }
    std::cout<<std::endl;
  }

  std::cout<<"Integrator settings: "<<std::endl;
  valueCheck(read_in_settings.dt,"dt");
  valueCheck(read_in_settings.t_f,"t_f");
  valueCheck(read_in_settings.bin_size,"bin_size");
  std::cout<<"Material properties: "<<std::endl;
  valueCheck(read_in_settings.k_n,"k_n");
  valueCheck(read_in_settings.k_t,"k_t");
  valueCheck(read_in_settings.eta_n,"eta_n");
  valueCheck(read_in_settings.eta_t,"eta_t");
  valueCheck(read_in_settings.m_mu,"mu");
  valueCheck(read_in_settings.m_rho,"rho");
  valueCheck(read_in_settings.freq,"freq");
  std::cout<<"out_forces: "<<read_in_settings.out_forces<<std::endl;


  return read_in_settings;

}
