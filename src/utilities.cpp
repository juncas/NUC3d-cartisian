#include "utilities.h"
#include <iostream>
#include <sstream>
#include <fstream>

void writeField_binary(std::ofstream &myFile, nuc3d::Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value=myField.getValue(i,j,k);
                myFile.write(reinterpret_cast<char *>(&value), sizeof(value));
            }
        }
    }
}


void readField_binary(std::ifstream &myFile, nuc3d::Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value;
                myFile.read(reinterpret_cast<char *>(&value), sizeof(value));
                myField.setValue(i,j,k,value);
            }
        }
    }
}


void writeField(std::ofstream &myFile, nuc3d::Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value=myField.getValue(i,j,k);
                myFile<<std::setprecision(12)<<value<<"\n";
            }
        }
    }
}


void readField(std::ifstream &myFile, nuc3d::Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value;
                if(!(myFile>>value))
                {
                    std::cout<<"Error:End of File at "
                    <<"i= "<<i<<", j="<<j<<", k= "<<k
                    <<std::endl;
                    exit(-1);
                }
                myField.setValue(i,j,k,value);
            }
        }
    }
}
