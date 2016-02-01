using msh
using PyPlot
using ODE

n = 10 #number of particles will be n*n
N = n*n
#set up mesh
m = msh.Mesh(n,n)
#msh.plotMesh(m,"b")
npull = 10 #number of elements being pulled at a constant speed
nhold = 10 #number of elements being held at their initial position

#Set up matrix of spring constants
global K
K = zeros(N,N)
for i = 1:N
  for j = i+1:N
    K[i,j] = 2^(-norm(m.coords[:,i] - m.coords[:,j],2)+1) #compute 2-norm because why not
    K[j,i] = K[i,j]
  end
end
#strengthen bond to the boundary nodes
for i = 1:nhold
  K[i,:] *= 10
  K[:,i] *= 10
end
for i = N-npull+1:N
  K[i,:] *= 10
  K[:,i] *= 10
end
#count the number of bonds each element has
global nbonds
nbonds = N*ones(N)

c1 = [4 2]  #choose speed in x and y direction
function f(t,x)
  println(t)
  global K
  #Create "forcing" term
  B = zeros(2*(N-npull-nhold),1)
  H = ones(size(K)[1],size(K)[2])
  xh  = [m.coords[1,1:nhold]' m.coords[2,1:nhold]'; hcat(x[1:N-npull-nhold], x[N-npull-nhold+1:2*(N-npull-nhold)]); (c1[1]*t+m.coords[1,N-npull+1:N])' (c1[2]*t+m.coords[2,N-npull+1:N])']'
  for i = 1:N
    for j = i+1:N
      if norm(xh[:,i]-xh[:,j]) > 10 
        H[i,j] = H[j,i] = 0
      end
      #make particles repulse each other/try to enforce order
      #TODO confuses the ode solver for some reason
      #if norm(xh[:,i] - xh[:,j]) < 0.1 && H[i,j] > 0
      #  H[i,j] = -1*H[i,j]
      #  H[j,i] = H[i,j]
      #end
    end 
  end
  K = K.*H
  #recalculate the number of bonds
  nbonds_t = sum((K .!= 0),2)
  global nbonds
  nbonds = [nbonds nbonds_t]

  #Create stiffness matrix
  A = zeros(N-npull-nhold,N-npull-nhold)
  for i = 1:N-npull-nhold
    A[i,i] = -sum(K[i+nhold,1:N]) #diagonal is special :)
    for j = 1:N-npull-nhold
      if(i!=j) A[i,j] = K[i+nhold,j+nhold] end
    end
  end
  for i = 1:N-npull-nhold
    B[i]     = sum((c1[1]*t + m.coords[1,N-npull+1:N]).*K[i+nhold,N-npull+1:N]) + sum(K[i+nhold, 1:nhold].*m.coords[1,1:nhold])
    B[i+N-npull-nhold] = sum((c1[2]*t + m.coords[2,N-npull+1:N]).*K[i+nhold,N-npull+1:N]) + sum(K[i+nhold, 1:nhold].*m.coords[2,1:nhold])
  end
    y = [A*x[1:N-npull-nhold]+B[1:N-npull-nhold]; A*x[N-npull-nhold+1:2*(N-npull-nhold)]+B[N-npull-nhold+1:2*(N-npull-nhold)]]'' 
  return [x[2*(N-npull-nhold)+1:end]; y]
end



maxT = 4 
y0 = [m.coords[1,nhold+1:N-npull]'; m.coords[2,nhold+1:N-npull]'; zeros(2*(N-npull-nhold))]
tout,yout  = ode45(f, y0,[0:0.1:maxT;])
ys = hcat(yout...)'
@printf "Number of timesteps: %d\n" size(tout)[1]

#Plot position of particles at time step
#count number of elements with broken bonds
broken = [nbonds .< 50]
cnt  = 1
minh = minimum(tout[2:end] - tout[1:end-1])
diff = maximum(tout[2:end] - tout[1:end-1])/minh 
@printf "Difference in time step size %d \n" diff
if diff > 5 i_step = floor(diff/5) end
for tstep = 1:i_step:size(tout)[1] - 1
  println(tstep)
  currenth = tout[tstep+1] - tout[tstep]
  for j =1:i_step:ceil(currenth/minh)
    broken_t = broken[:,tstep]
    yh  = [m.coords[1,1:nhold]' m.coords[2,1:nhold]';
           vcat(ys[tstep, 1:N-npull-nhold], ys[tstep, N-npull-nhold+1:2*(N-npull-nhold)])';
           (c1[1]*tout[tstep]+m.coords[1,N-npull+1:N])' (c1[2]*tout[tstep]+m.coords[2,N-npull+1:N])']
    plot(yh[ broken_t,1], yh[ broken_t,2], marker="o", markersize=10, color = "r", linestyle = "None")
    #plot(yh[~broken_t,1], yh[~broken_t,2], marker="o", markersize=10, color = "b", linestyle = "None")
    savefig(string("Plots/plot", cnt))
    clf()
    cnt = cnt+1
  end
end

