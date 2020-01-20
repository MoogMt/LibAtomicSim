include("xyz.jl")

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/9.8/3000K/"
atoms = filexyz.read(string(folder,"first.xyz"),1)[1]

vol=9.8

targetVol=9.7

cell=cell_mod.Cell_param(targetVol,targetVol,targetVol)

for atom1=1:size(atoms.names)[1]
    for atom2=atom1+1:size(atoms.names)[1]
        distanceatm = cell_mod.distance(atoms,cell,atom1,atom2)
        if distanceatm < 1.0
            print("atom1: ",atom1," atom2: ",atom2," distance:",distanceatm,"\n")
        end
    end
end

Cfile=open(string("/home/moogmt/smallC_",targetVol,".cpmd"),"w")
Ofile=open(string("/home/moogmt/smallO_",targetVol,".cpmd"),"w")
for atom=1:size(atoms.names)[1]
    if atoms.names[atom] == "C"
        for i=1:3
            write(Cfile,string(atoms.positions[atom,i]," "))
        end
        write(Cfile,"\n")
    else
        for i=1:3
            write(Ofile,string(atoms.positions[atom,i]," "))
        end
        write(Ofile,"\n")
    end
end
close(Cfile)
close(Ofile)
