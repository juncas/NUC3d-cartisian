#ifndef utilities_h
#define utilities_h
#include <memory>
#include "field.h"


void writeField_binary(std::ofstream &myFile, nuc3d::Field &myField);
void readField_binary(std::ifstream &myFile, nuc3d::Field &myField);
void writeField(std::ofstream &myFile, nuc3d::Field &myField);
void readField(std::ifstream &myFile, nuc3d::Field &myField);


#endif
