#include <iostream>
#include <fstream>
#include <vector>
//
// Created by maytee on 9/18/19.
//

#ifndef DEM_OUTPUTWRITER_H
#define DEM_OUTPUTWRITER_H

void writeRegionVTU(std::string file_name, std::vector<double>& grid_min, std::vector<double>& grid_max);
void writeGrainsVTU(std::string file_name, std::vector<double>& q, std::vector<double>& v, std::vector<double>& r,
                    std::vector<int>& unique_id);

std::string generateOutputFileStr(std::string prefix, int freq, double t_f, double dt, int& output_num);



#endif //DEM_OUTPUTWRITER_H
