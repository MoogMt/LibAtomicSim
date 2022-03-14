// Guard
#ifndef frame_h
#define frame_h

// Frame for a given structure
class Frame 
{
  protected:
    // Class variables
    //-------------------------
    int number;
    std::string* names;
    double* positions; 
    double* cell_matrix;
    //-------------------------
  
  public:
    // Constructors
    Frame();
    Frame( int new_number );
    Frame( int new_number, std::string* new_names, double* new_positions );

    // Getters
    int getNumberAtoms();
    std::string* getNames();
    double* getPositions();
    double* getCellMatrix();

    // Setters
    // -> Number of atoms
    void setNumberAtoms( int new_number );
    // -> Names
    void setNames( std::string* new_names );
    void setName( std::string new_name, int index_atom );
    void setNames( double* new_names, int* index_atoms, int number_new );
    // -> Positions
    void initPosition( double* new_positions, int new_number );
    void setPositions( double* new_positions );
    void setPosition( double* new_positions, int atom_index );
    void setPositions( double* new_positions, int* atom_indexs, int number_new );
  
    // Modifyers
    void addAtoms( int number_to_add, std::string* new_names, double* new_positions );

    // Checkers
    void printAtoms();
};

#endif
