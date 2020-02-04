module pimaim

using LinearAlgebra

using conversion
using atom_mod
using cell_mod
using filexyz
using pdb
using periodicTable

function writeRestart( path_file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_matrix }
    nb_atoms=atom_mod.getNbAtoms(atoms)
    handle_out=open( path_file , "w" )
    write( handle_out, string("T\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    for atom=1:nb_atoms
        for i=1:3
            write( handle_out, string( atoms.positions[atom,i]*conversion.ang2Bohr," " ) )
        end
        write( handle_out, string("\n") )
    end
    matrix=cell.matrix
    norms_v=zeros(Real,3)
    for i=1:3
        norms_v[i] = LinearAlgebra.norm(matrix[i,:])
        matrix[:,i] = matrix[:,i]./norms_v[i]
    end
    for i=1:3
        for j=1:3
            if matrix[i,j] < 10^(-5)
                matrix[i,j] = 0
            end
            write( handle_out, string( matrix[i,j], " ") )
        end
        write( handle_out, string("\n") )
    end
    for i=1:3
        write(handle_out,string(norms_v[i]*conversion.ang2Bohr,"\n"))
    end
    close(handle_out)
end

function writeCrystalCell( path_file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_matrix }

    norms_v=zeros( Real, 3 )
    matrix = copy( cell.matrix )
    species = atom_mod.getSpecies( atoms )
    for i=1:3
        norms_v[i] = LinearAlgebra.norm(matrix[i,:])
        matrix[:,i] = matrix[:,i]./norms_v[i]
    end

    handle_in = open( path_file, "w" )
    for i=1:3
        for j=1:3
            if matrix[j,i] < 0.001
                write( handle_in, string( 0 , " " ) )
            else
                write( handle_in, string( matrix[j,i], " " ) )
            end
        end
        write( handle_in, string("\n") )
    end
    for i=1:3
        write( handle_in, string( norms_v[i], "\n")  )
    end
    for i=1:3
        write( handle_in, "1\n" )
    end
    for i=1:size(species)[1]
        write( handle_in, string( periodicTable.names2Z( species[i] ), "\n" ) )
        write( handle_in, string( species[i], "_quartz.mat\n" ) )
    end
    for i=1:3
        write( handle_in, string( norms_v[i], "\n")  )
    end
    close(handle_in)
    return true
end

end
