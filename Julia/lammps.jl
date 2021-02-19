module lammps

# Loading necessary modules from LibAtomSim
using conversion
using atom_mod

# Description
# Handles interface with LAMMPS: Reads input and output files, writes input file,
# sets up simulations based on user parameters, etc...

export writeLammpsInput

# Writting from a collection of molecules
#-------------------------------------------------------------------------------
function writeLammpsInput( molecules::Vector{T1}, file_path::T2 ) where { T1 <: atom_mod.AtomMolList, T2 <: AbstractString }
    # Argument
    # - molecules: AtomMolList containing information about atoms and molecules
    # - file_path: path to the output file
    # Output
    # - Bool: whether the writting was successful

    # Get number of types
    nb_types = size( getNames( molecules ) )[1]

    # Get number of atoms
    nb_atoms = getAtomsNb( molecules )

    # Opening output file
    handle_out = open( file_path, "w+" )

    # Writes first line of file
    write( handle_out, "LAMMPS data file, timestep = 0\n" )

    # Writes empty line
    write( handle_out, "\n" )

    # write number of atoms
    write( handle_out, string( nb_atoms, " atoms\n" ) )

    # Write number of different atomic types
    write( handle_out, string( nb_types, " atom types\n" ) )

    # Write empty line
    write( handle_out,"\n")

    # Get X,Y,Z positions
    X = getX( molecules )
    Y = getY( molecules )
    Z = getZ( molecules )

    # Writes minimum and maximum values in each direction
    write( handle_out, string( getMin(X), " ", getMax(X), " xlo xhi\n" ) )
    write( handle_out, string( getMin(Y), " ", getMax(Y), " ylo yhi\n" ) )
    write( handle_out, string( getMin(Z), " ", getMax(Z), " zlo zhi\n" ) )

    # Empty cartesian positions (X,Y,Z)
    empty!( X )
    empty!( Y )
    empty!( Z )

    # Writes empty line again
    write(out,"\n")

    # Write "Masses", start masses section
    write(out,"Masses\n")

    # Writes empty line again
    write(out,"\n")

    # Get masses of atoms
    masses=getMasses(molecules)

    # Get labels of atoms
    labels=getLabels(molecules)

    # Loop over types
    for type=1:nb_types
        # Write the masses for each atoms label
        write( handle_out, string( labels[type], " ", masses[type], "\n" ) )
    end

    # Writes empty line again
    write(out,"\n")

    # Write "ATOMS", start Atoms section
    write(out,"Atoms\n")

    # Writes empty line again
    write(out,"\n")

    # Start atom counter to 1
    atom_nb=1

    # Get number of molecules
    nb_molecules = size(molecules)[1]

    # Loop over molecules
    for molecule=1:nb_molecules
        # Get the atoms from a given molecule
        atoms = getAtoms( molecules[molecule] )

        # Loop over atoms
        for atom=1:size(atoms)[1]
            # Writes atom index
            write( handle_out, string( atom_nb, " " ) )
            # Writes molecule index
            write( handle_out, string( molecule, " " ) )
            # Write atom label
            write( handle_out, string( getLabel( atoms[j] ), " " ) )
            # Write atom charge
            write( handle_out, string( getCharge (atoms[j] )," " ) )
            # Write atomic positions
            write( handle_out, string( getX( atoms[j] ), " ", getY( atoms[j] ), " ", getZ( atoms[j] ) ) )
            # Writes end of line
            write( handle_out, "\n")
            # Increments atom number
            atom_nb = atom_nb + 1
        end
    end

    # Closing output file
    close(handle_out)

    # Returns true if writting was successful
    return true
end
#-------------------------------------------------------------------------------

end
