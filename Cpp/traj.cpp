#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "utils.h"
#include "complex_numbers.h"
#include "linear_algebra.h"
#include "frame.h"
#include "traj.h"

using namespace std;

Traj::Traj()
{
    number = 0;
    names_code=(int*) malloc(0);
    positions_cart  = zeroMatrix(3,0); 
    positions_red   = zeroMatrix(3,0);
    cell_matrix     = initCellMatrix();
    inv_cell_matrix = initCellMatrix();
}
