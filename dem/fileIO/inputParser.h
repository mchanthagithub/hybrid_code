//
// Created by maytee on 11/14/19.
//

#include <string>
#include <iostream>
#include <fstream>
#include "../InputSettings.h"

#ifndef DEM_INPUTPARSER_H
#define DEM_INPUTPARSER_H

template<typename T> T getValueFromEqualityString(const char* value_char, std::string line);


InputSettings readInputFile(std::string file_name);


#endif //DEM_INPUTPARSER_H
