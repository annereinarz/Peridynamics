module msh

using PyPlot

export plotMesh

struct Mesh
    coords

    function Mesh(n1,n2)
      @assert (n1 > 0 && n2 > 0) "Mesh arguments must be positive"
      coords = createCoords(n1,n2)
      new(coords)
    end
end

  function createCoords(n1,n2)

  #Fill coords array
  coords = zeros(2,n1*n2);
  cnt = 1
  for i = 1:n2
    for j = 1:n1
      coords[1:2,cnt] = [j-1,i-1]
      cnt = cnt+1
     end
  end
  r = 1.0/(4*n1)*rand(size(coords)[1], size(coords)[2])
  coords  = coords + r  #"random" pertubation of grid
  return coords
  end


  function plotMesh(m::Mesh,c::String)

    #Plot coordinates
    for i = 1:size(m.coords)[2]
        plot(m.coords[1,i],m.coords[2,i],marker="o",color =c, markersize = 10)
    end

  end
end
