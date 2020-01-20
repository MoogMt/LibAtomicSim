#ifndef GEOM_H
#define GEO_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

struct Point
{
  double x, double y, double z;
}

struct Points
{
  std::vector<double> x, y, z;
};

#endif
