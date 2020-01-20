# getColvar

function getColvar(file::AbstractString, stepinit::Int, stepend::Int, nb_cvs::Int )
  #------------------------------------------------------------------------------
  # in:
  #  -> file (string) : colvar path
  #  -> stepinit (int) : begin step
  #  -> stepend (int): end step
  #  -> nb_cvs (int) : number of collective variables
  #----------------------------------------------------
  # out:
  # -> timestep:
  #----------------------------------------------------------
  timestep=[]; cvs=zeros(Float64,1,nb_cvs); bias=[]
  open(file,"r") do f
    for line in eachline(f)
      if ( ! contains(line,"#") )
        step=parse(Float64,split(line)[1])
        if ( step >= stepinit && step <= stepend )
          timestep=[timestep;step]
          temp=zeros(Float64,1,nb_cvs)
          for i=1:nb_cvs
            temp[1,i]=parse(Float64,split(line)[1+i])
          end
          cvs=[cvs;temp]
          push!(bias,parse(Float64,split(line)[nb_cvs+2]))
        end
      end
    end
  end
  cvs=cvs[2:size(cvs)[1],:]
  return timestep,cvs,bias
end
