#include <iostream>
#include <fstream>
#include <vector>
#include "../Region.h"
//
// Created by maytee on 9/18/19.
//

#ifndef DEM_OUTPUTWRITER_H
#define DEM_OUTPUTWRITER_H

void writeRegionVTU(std::string file_name, std::vector<double>& grid_min, std::vector<double>& grid_max);
void writeGrainsVTU(std::string file_name, std::vector<double>& q, std::vector<double>& v, std::vector<double>& r,
                    std::vector<int>& unique_id, std::vector<double>& f);
void writeForcesVTU(std::string file_name, std::map<std::pair<int,int>,ContactGrainGrain> grain_grain_contacts,
                    std::map<std::pair<int,int>,ContactGrainWall> grain_wall_contacts);

std::string generateOutputFileStr(std::string prefix, int freq, double t_f, double dt, int& output_num);



#endif //DEM_OUTPUTWRITER_H
