// Guard
#ifndef traj_h
#define traj_h

// Object for trajectories
class Traj
{
  protected:
    // Class variables
    //-------------------------
    int number;
    int* names_code;
    double* positions_cart; 
    double* positions_red;
    double* cell_matrix;
    double* inv_cell_matrix;
    //-------------------------

  public:
    // Constructors
    Traj();
    //Traj( int nb_atoms );
    //Traj( int nb_atoms, int, nb_steps );
    //Traj( int new_number, std::string* new_names, double* new_positions );
};


#endif