module plumed

# Rational Switching Functions
#------------------------------------------------------------------------------------------------
# Forms:
#               1 - [(x-d_0)/r_0]^n
# -> swf(r) = -----------------------
#               1 - [(x-d_0)/r_0]^m
#
mutable struct switchingFullRational
    #------------------------------------------
    # Parameters of the structure
    #------------------------------------------
    r_0::Real # Reference distance in nm
    d_0::Real # Offset distance in nm
    n::Int    # Integer
    m::Int    #
    #------------------------------------------
    # Constructors
    #------------------------------------------
    function switchingFullRational()
        new( 0.1, 0.0, 1, 2 );
    end
    function switchingFullRational( r_0::T1, n::T2 ) where { T1 <: Real, T2 <: Int }
        new( r_0, 0.0, n, 2*n )
    end
    function switchingFullRational( r_0::T1, n::T2, m::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
        if n >= m
            print("ERROR n >= m when n < m is expected.\n")
            print("n = ",n," m=",m,"\n")
            return false
        end
        new( r_0, 0.0, n, m )
    end
    function switchingFullRational( r_0::T1, d_0::T2, n::T3, m::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Int }
        if n >= m
            print("ERROR n >= m when n < m is expected.\n")
            print("n = ",n," m=",m,"\n")
            return false
        end
        new( r_0, d_0, n, m )
    end
    #------------------------------------------
end
function evaluate( swf::T1, r::T2 ) where { T1 <: switchingFullRational, T2 <: Real }
    return (1-((r-swf.d_0)/swf.r0)^swf.n)/(1-((r-swf.d_0)/swf.r0)^swf.m)
end
function evaluate( swf::T1, r::Vector{T2} ) where { T1 <: switchingFullRational, T2 <: Real }
    results=zeros( size(r)[1] )
    for i=1:size(r)[1]
        results[i] = evaluate( swf, r[i] )
    end
    return results
end
function evaluate( swf::T1, r::Array{T2,2} ) where { T1 <: switchingFullRational, T2 <: Real }
    results=zeros( size(r)[1], size(r)[2] )
    for i=1:size(r)[1]
        for j=1:size(r)[2]
            results[i,j] = evaluate( swf, r[i,j] )
        end
    end
    return results
end
function writePlumedInput( handle_out::T1, swf::T2 ) where { T1 <: IO, T2 <: switchingFullRational }
    Base.write( handle_out, string( "RATIONAL ") )
    Base.write( handle_out, string( "R_0=",swf.r_0," ") )
    Base.write( handle_out, string( "D_0=",swf.d_0," ") )
    Base.write( handle_out, string( "MM=",swf.r_m," ") )
    Base.write( handle_out, string( "NN=",swf.n," ") )
    return true
end
# Same as above, but with m=2*n so that the evaluation is faster with the trick:
#                     1
# -> swf(r) = -----------------------
#               1 + [(x-d_0)/r_0]^n
mutable struct switchingRational
    #------------------------------------------
    # Parameters of the structure
    #------------------------------------------
    r_0::Real
    d_0::Real
    n::Int
    #------------------------------------------
    # Constructors
    #------------------------------------------
    function switchingRational()
        new( 0.1, 0.0, 1, 2 );
    end
    function switchingRational( r_0::T1 ) where { T1 <: Real, T2 <: Int }
        new( r_0, 0.0, 1 )
    end
    function switchingRational( r_0::T1, n::T2 ) where { T1 <: Real, T2 <: Int }
        new( r_0, 0.0, n )
    end
    function switchingRational( r_0::T1, d_0::T2 ) where { T1 <: Real, T2 <: Real }
        new( r_0, d_0, 1 )
    end
    function switchingRational( r_0::T1, d_0::T2, n::T2 ) where { T1 <: Real, T2 <: Real, T3 <: Int }
        new( r_0, d_0, n )
    end
    #------------------------------------------
end
function evaluate( swf::T1, r::T2 ) where { T1 <: switchingRational, T2 <: Real }
    return 1/(1+((r-swf.d_0)/r_0)^(swf.n))
end
function evaluate( swf::T1, r::Vector{T2} ) where { T1 <: switchingRational, T2 <: Real }
    results=zeros( size(r)[1] )
    for i=1:size(r)[1]
        results[i] = evaluate( swf, r[i] )
    end
    return results
end
function evaluate( swf::T1, r::Array{T2,2} ) where { T1 <: switchingRational, T2 <: Real }
    results=zeros( size(r)[1], size(r)[2] )
    for i=1:size(r)[1]
        for j=1:size(r)[2]
            results[i,j] = evaluate( swf, r[i,j] )
        end
    end
    return results
end
function writePlumedInput( handle_out::T1, swf::T2 ) where { T1 <: IO, T2 <: switchingRational }
    Base.write( handle_out, string( "RATIONAL ") )
    Base.write( handle_out, string( "R_0=",swf.r_0," ") )
    Base.write( handle_out, string( "D_0=",swf.d_0," ") )
    Base.write( handle_out, string( "MM=",Int(swf.n*2)," ") )
    Base.write( handle_out, string( "NN=",swf.n," ") )
    return true
end
# Same as rational, but with d_0=0.0
mutable struct switchingSimpleRational
    #------------------------------------------
    # Parameters of the structure
    #------------------------------------------
    r_0::Real
    n::Int
    #------------------------------------------
    # Constructors
    #------------------------------------------
    function switchingSimpleRational()
        new( 0.1, 1 );
    end
    function switchingSimpleRational( r_0::T1 ) where { T1 <: Real, T2 <: Int }
        new( r_0, 1 )
    end
    function switchingSimpleRational( r_0::T1, n::T2 ) where { T1 <: Real, T2 <: Int }
        new( r_0, n )
    end
    #------------------------------------------
end
function evaluate( swf::T1, r::T2 ) where { T1 <: switchingSimpleRational, T2 <: Real }
    return 1/(1+(r/swf.r_0)^(swf.n))
end
function evaluate( swf::T1, r::Vector{T2} ) where { T1 <: switchingSimpleRational, T2 <: Real }
    results=zeros( size(r)[1] )
    for i=1:size(r)[1]
        results[i] = evaluate( swf, r[i] )
    end
    return results
end
function evaluate( swf::T1, r::Array{T2,2} ) where { T1 <: switchingSimpleRational, T2 <: Real }
    results=zeros( size(r)[1], size(r)[2] )
    for i=1:size(r)[1]
        for j=1:size(r)[2]
            results[i,j] = evaluate( swf, r[i,j] )
        end
    end
    return results
end
function writePlumedInput( handle_out::T1, swf::T2 ) where { T1 <: IO, T2 <: switchingSimpleRational }
    Base.write( handle_out, string( "RATIONAL ") )
    Base.write( handle_out, string( "R_0=",swf.r_0," ") )
    Base.write( handle_out, string( "D_0=",0.0," ") )
    Base.write( handle_out, string( "MM=",Int(swf.n*2)," ") )
    Base.write( handle_out, string( "NN=",swf.n," ") )
    return true
end
#------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
switchingFunction=Union{ switchingFullRational, switchingRational, switchingSimpleRational }
function writeSwitchFunction( handle_out::T1, number::T2, swf::T3 ) where { T1<:IO, T2 <: Int, T3 <: switchingFunction }
    Base.write( handle_out, string( "SWITCH", number,"={ " ) )
    writePlumedInput( handle_out, swf )
    Base.write( handle_out, string( " }\n" ) )
    return true
end
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
mutable struct pathZS
    #----------------------------------------------------------
    label::AbstractString
    args::Vector{AbstractString}
    lambda::Real
    #----------------------------------------------------------
    function pathZS( label::T1, args::Vector{T2} ) where { T1 <: AbstractString, T2 <: AbstractString }
        new( label, args, 1.0 )
    end
    function pathZS( label::T1, args::Vector{T2}, lambda::T3 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real }
        new( label, args, lambda )
    end
    #----------------------------------------------------------
end
function writeActDispatch( handle_out, act2 ) where { T1 <: IO, T2 <: pathZS }
    Base.write( handle_out, string(pathZS.label,": ") )
    Base.write( handle_out, string("FUNCPATHMSD ") )
    #------------------------------------------------
    Base.write( handle_out, string("ARG= ") )
    nb_args=size(args)[1]
    for arg=1:nb_args
        Base.write( handle_out, string( pathZS.args[arg] ) )
        if arg < nb_args
            Base.write( handle_out, string(",") )
        end
    end
    Base.write( handle_out, string( " " ) )
    #------------------------------------------------
    Base.write( handle_out, string( "LAMBDA=", pathZS.lambda ) )
    return true
end
plumedAct=Union{ pathZS }
#------------------------------------------------------------------------------------------------

# INPUT handling
#------------------------------------------------------------------------------------------------
mutable struct inputPIV
    #------------------------------------------
    # Parameters of the structure
    #------------------------------------------
    label::AbstractString
    ref_files::Vector{AbstractString}
    atom_names::Vector{AbstractString}
    switching_fct::Vector{switchingFunction}
    cut_off::Real
    skin::Real
    stride::Int
    s_factors::Vector{Real}
    volume::Bool
    onlydirect::Bool
    #------------------------------------------
    # Constructors
    #------------------------------------------
    function inputPIV( label::T1,
        ref_files::Vector{T2},
        atom_names::Vector{T3},
        switching_fct::Vector{T4},
        cut_off::T5,
        skin::T6,
        stride::T7,
        s_factors::Vector{T8},
        volume::T9,
        onlydirect::T10 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: switchingFunction, T5 <: Real, T6 <: Real, T7 <: Int , T8 <: Real, T9 <: Bool, T10 <: Bool }
        new( label, ref_files, atom_names, switching_fct, cut_off, skin, stride, s_factors, volume, onlydirect )
    end
end
function writeInputPIV( handle_out::T1, input_::T2 ) where { T1 <: IO, T2 <: inputPIV }
    Base.write( handle_out, string( "PIV ...", "\n" ) )
    # Label for the descriptor
    #----------------------------------------------------------
    Base.write( handle_out, string( "LABEL=", input_.label , "\n") )
    # Choosing whether or not to use volume rescale
    #----------------------------------------------------------
    if input_.volume
        Base.write( handle_out, string( "VOLUME", "\n") )
    end
    # Path tot the references
    #----------------------------------------------------------
    nb_ref = size( input_.ref_files )[1]
    for ref=1:nb_ref
        Base.write( handle_out, string( "REF_FILE",ref,"=", input_.ref_files[ref] , "\n") )
    end
    # Target atoms species to compute PIV
    #----------------------------------------------------------
    nb_atoms = size( input_.atom_names )[1]
    Base.write( handle_out, string( "ATOMTYPES=",) )
    for atom=1:nb_atoms
        Base.write( handle_out, string( input_.atom_names[atom] ) )
        if atom < nb_atoms
            Base.write( handle_out, string(","))
        end
    end
    Base.write( handle_out, string("\n") )
    # S factors, allows to modulate the importance of the block
    #----------------------------------------------------------
    nb_block=size(input_.s_factors)[1]
    Base.write( handle_out, string( "SFACTOR=",) )
    for iblock=1:nb_block
        Base.write( handle_out, string( input_.s_factors[iblock] ) )
        if iblock < nb_block
            Base.write( handle_out, string(","))
        end
    end
    Base.write( handle_out, string("\n") )
    # Switching functions, one per block
    #----------------------------------------------------------
    for iblock=1:nb_block
        writeSwitchFunction( handle_out, iblock, input_.switching_fct[iblock] )
    end
    # Cut Off (nm)
    #----------------------------------------------------------
    Base.write( handle_out, string( "NL_CUTOFF=", input_.cut_off , "\n") )
    # stride ?
    #----------------------------------------------------------
    Base.write( handle_out, string( "NL_STRIDE=", input_.stride , "\n") )
    # Skin (nm)
    #----------------------------------------------------------
    Base.write( handle_out, string( "NL_SKIN=", input_.skin , "\n") )
    # Not counting cross blocks (O-H distances for water for example)
    #----------------------------------------------------------
    if input_.onlydirect
        Base.write( handle_out, string( "ONLYDIRECT", "\n") )
    end
    #----------------------------------------------------------
    Base.write( handle_out, string( "... PIV", "\n" ) )
end
function writeAct( handle_out::T1, act::T2 ) where { T1 <: IO, T2 <: plumedAct }
    return writeActDispatch( handle_out, act2 )
end
function writeOutputInstruction( handle_out::T1 , args::Vector{T2}, stride::T3, file_out_name::T4="COLVAR", format::T5="%15.6f" ) where { T1 <: IO, T2 <: AbstractString, T3 <: Int, T4 <: AbstractString, T5 <: AbstractString }
    Base.write( handle_out, string("PRINT ") )
    #-----------------------------------------
    Base.write( handle_out, string("ARG=") )
    nb_arg=size(args)
    for arg=1:nb_arg
        Base.write( handle_out, string( args[i] ) )
        if arg < nb_arg
            Base.write( handle_out, string(",") )
        end
    end
    #-----------------------------------------
    Base.write( handle_out, string( "STRIDE=", stride, " " ) )
    #-----------------------------------------
    Base.write( handle_out, string( "FILE=", file_out_name, " " ) )
    #-----------------------------------------
    Base.write( handle_out, string( "FMT=", format, " " ) )
    #-----------------------------------------
    return true
end
#------------------------------------------------------------------------------------------------



# COLVAR handling
#------------------------------------------------------------------------------------------------
function getNbColumnColvar( file_path::T1 ) where { T1 <: AbstractString }
    handle_in = open( file_path, "w" )
    readline( handle_in)
    nb_col=size( split( readline(handle_in) ) )[1]
    close( handle_in )
    return nb_col
end
function readColvar( file_path::T1 ) where { T1 <: AbstractString }
    # Getting nb of lines
    nb_lines = getNbLines( file_path )
    if nb_lines == false
        return false
    end
    # Columns
    nb_col = getNbColumnColvar( file_path )
    if nb_col < 2
        print("Problem with colvar file at: ",file_path,"\n")
        print("Number of column to low: nb_col=",nb_col,"\n")
        print("We expect at least 2 colums (timestep + variable)\n")
        return false
    end
    time_ = zeros( nb_lines )
    data_ = zeros( nb_lines, nb_col-1 )
    handle_in = open( file_path )
    readline( handle_in )
    for line=1:nb_lines
        keys = split( readline(handle_in) )
        time_ = parse(Float64, keys[1] )
        for col_dat=1:nb_col-1
            data[line,col_dat] = parse(Float64, keys[col_dat+1])
        end
    end
    close( handle_in )
    return time,data
end
#------------------------------------------------------------------------------------------------

end
