using Plots
using LinearAlgebra

# setting parameters
Reynolds_number = 500 ;
max_time_step = 50 ;
length_of_x = 1 ;
length_of_y = 1 ;
number_of_grids_in_x = 80 ;
number_of_grids_in_y = 80 ;
Δt = 0.01 ;

width_of_grid_in_x = length_of_x / number_of_grids_in_x ;   # dx
width_of_grid_in_y = length_of_y / number_of_grids_in_y ;   # dy

# define initial physical quantity
u = zeros(number_of_grids_in_x+1, number_of_grids_in_y+2) ; # velocity in x direction (u)
v = zeros(number_of_grids_in_x+2, number_of_grids_in_y+1) ; # velocity in y direction (v)
p = zeros(number_of_grids_in_x, number_of_grids_in_y) ; # pressure

function solvePoissonEquation_direct(RHS, Nx, Ny, dx, dy)
    Abig = zeros(Nx*Ny, Nx*Ny) ;
    Ax = zeros(Nx, Nx) ;
    for row = 1:Nx
        if row < Nx
            Ax[row, row+1] = 1 ;
        end
        Ax[row, row] = -2 ;
        if row > 1
            Ax[row, row-1] = 1 ;
        end
    end
    Ax[1, 1] = -1 ;
    Ax[Nx, Nx] = -1 ;
    
    dd = Matrix(1.0I, Nx, Nx);
    for jj=1:Ny
        if jj==1 || jj==Ny
            Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-1):Nx*jj] .= Ax ./ dx.^2 .- dd ./ dy.^2 ;
        else
            Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-1):Nx*jj] .= Ax ./ dx.^2 .- 2 .* dd ./ dy.^2 ;
        end
        if jj!=1
            Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-2):Nx*(jj-1)] .= dd ./ dy.^2;
        end
        if jj!=Ny
            Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-0):Nx*(jj+1)] .= dd ./ dy.^2;
        end
    end

    Abig[1,:] .= 0 ;
    Abig[1,1] = 1 ;
    RHS[1] = 0 ;

    u = Abig\RHS ;
    u = reshape(u, Nx, Ny) ;

    return u
end

for K = 1:max_time_step
    bctop = 1; # 境界上部の速度 u
    u[:, 1] .= - u[:, 2];
    v[:, 1] .= 0;             # bottom
    u[:,size(u)[2]] .= 2 .* bctop .- u[:, size(u)[2]-1];
    v[:, size(v)[2]] .= 0;  # top
    u[1, :] .= 0;
    v[1, :] .= -v[2, :];             # left
    u[size(u)[1], :] .= 0;
    v[size(v)[1], :] .= -v[size(v)[1]-1, :];    # right

    ∇∇ux = (u[1:size(u)[1]-2, 2:size(u)[2]-1] .- 2 .* u[2:size(u)[1]-1, 2:size(u)[2]-1] .+ u[3:size(u)[1], 2:size(u)[2]-1]) ./ width_of_grid_in_x.^2 ; 
    ∇∇uy = (u[2:size(u)[1]-1, 1:size(u)[2]-2] .- 2 .* u[2:size(u)[1]-1, 2:size(u)[2]-1] .+ u[2:size(u)[1]-1, 3:size(u)[2]]) ./ width_of_grid_in_y.^2 ;
    ∇∇vx = (v[1:size(v)[1]-2, 2:size(v)[2]-1] .- 2 .* v[2:size(v)[1]-1, 2:size(v)[2]-1] .+ v[3:size(v)[1], 2:size(v)[2]-1]) ./ width_of_grid_in_x.^2 ;
    ∇∇vy = (v[2:size(v)[1]-1, 1:size(v)[2]-2] .- 2 .* v[2:size(v)[1]-1, 2:size(v)[2]-1] .+ v[2:size(v)[1]-1, 3:size(v)[2]]) ./ width_of_grid_in_y.^2 ;

    # 1. interpolate velocity at cell center/cell cornder
    u_in_cell_center = (u[1:size(u)[1]-1, 2:size(u)[2]-1] .+ u[2:size(u)[1], 2:size(u)[2]-1]) ./ 2 ;
    u_in_cell_corner = (u[:, 1:size(u)[2]-1] .+ u[:, 2:size(u)[2]]) ./ 2 ;
    v_in_cell_center = (v[2:size(v)[1]-1, 1:size(v)[2]-1] .+ v[2:size(v)[1]-1, 2:size(v)[2]]) ./ 2 ;
    v_in_cell_corner = (v[1:size(v)[1]-1, :] .+ v[2:size(v)[1], :]) ./ 2 ;
    
    # 2. multiply
    uuce = u_in_cell_center .* u_in_cell_center ;
    uvco = u_in_cell_corner .* v_in_cell_corner ;
    vvce = v_in_cell_center .* v_in_cell_center ;
    
    # 3-1. get derivative for u
    Nu = (uuce[2:size(uuce)[1], :] .- uuce[1:size(uuce)[1]-1, :]) ./ width_of_grid_in_x ;
    Nu .+= (uvco[2:size(uvco)[1]-1, 2:size(uvco)[2]] .- uvco[2:size(uvco)[1]-1, 1:size(uvco)[2]-1]) ./ width_of_grid_in_y ;

    # 3-2. get derivative for v
    Nv = (vvce[:, 2:size(vvce)[2]] .- vvce[:, 1:size(vvce)[2]-1]) ./ width_of_grid_in_y ;
    Nv .+= (uvco[2:size(uvco)[1], 2:size(uvco)[2]-1] .- uvco[1:size(uvco)[1]-1, 2:size(uvco)[2]-1]) ./ width_of_grid_in_x ;

    # Get intermidiate velocity
    u[2:size(u)[1]-1, 2:size(u)[2]-1] .= u[2:size(u)[1]-1, 2:size(u)[2]-1] .+ Δt .* (-Nu .+ (∇∇ux .+ ∇∇uy) ./ Reynolds_number) ;
    v[2:size(v)[1]-1, 2:size(v)[2]-1] .= v[2:size(v)[1]-1, 2:size(v)[2]-1] .+ Δt .* (-Nv .+ (∇∇vx .+ ∇∇vy) ./ Reynolds_number) ;

    # velocity correction
    b = (u[2:size(u)[1], 2:size(u)[2]-1] .- u[1:size(u)[1]-1, 2:size(u)[2]-1]) ./ width_of_grid_in_x .+ (v[2:size(v)[1]-1, 2:size(v)[2]] .- v[2:size(v)[1]-1, 1:size(v)[2]-1]) ./ width_of_grid_in_y ;
    RHS = reshape(b, number_of_grids_in_x * number_of_grids_in_y, 1) ;
    
    p = solvePoissonEquation_direct(RHS, number_of_grids_in_x, number_of_grids_in_y, width_of_grid_in_x, width_of_grid_in_y)
    
    # 圧力勾配を仮の速度場から引いて、新しい速度場
    u[2:size(u)[1]-1, 2:size(u)[2]-1] .= u[2:size(u)[1]-1,2:size(u)[2]-1] .- (p[2:size(p)[1],:] .- p[1:size(p)[2]-1,:]) ./ width_of_grid_in_x ;
    v[2:size(v)[1]-1, 2:size(v)[2]-1] .= v[2:size(v)[1]-1,2:size(v)[2]-1] .- (p[:,2:size(p)[2]] .- p[:,1:size(p)[2]-1]) ./ width_of_grid_in_y ;
end

# output
u_in_cell_center = (u[1:size(u)[1]-1, 2:size(u)[2]-1] .+ u[2:size(u)[1], 2:size(u)[2]-1]) ./ 2 ;
v_in_cell_center = (v[2:size(v)[1]-1, 1:size(v)[2]-1] .+ v[2:size(v)[1]-1, 2:size(v)[2]]) ./ 2 ;

velocity_abs = sqrt.(u_in_cell_center.^2 .+ v_in_cell_center.^2) ;
contour(velocity_abs, fill = true, c = cgrad(:grays, rev = true)) ;