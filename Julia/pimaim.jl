module pimaim

using LinearAlgebra

using conversion
using atom_mod
using cell_mod
using filexyz
using pdb

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
        norms_v[i] = LinearAlgebra.norm(matrix[:,i])
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

end
