module periodicTable

# Description
# Functions that helps transform Z number to element name and vice
# versa

# Exporing functions
export z2Names
export names2Z

zNames=[ "H",                                                                                                                                                                                     "He",
        "Li", "Be",                                                                                                                                                  "B",  "C",  "N",  "O",  "F", "Ne",
        "Na", "Mg",                                                                                                                                                 "Al", "Si",  "P",  "S", "Cl", "Ar",
         "K", "Ca",                                                                                     "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr",                                                                                      "Y", "Zr", "Nb", "Mo", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sb", "Te",  "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

function z2Names( z::T1 ) where { T1 <: Int }
    return zNames[z]
end

function names2Z( name::T1 ) where { T1 <: AbstractString }
    for i=1:size(zNames)[1]
        if name == zNames[i]
            return i
        end
    end
    print("No such element in the periodic table\n")
    return false
end

end
