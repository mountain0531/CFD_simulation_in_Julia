using Plots
using LinearAlgebra
using BenchmarkTools

module Parameters
    struct Para
        Reynolds_number::Float64
        max_time_step::Int64
        length_of_x::Float64
        length_of_y::Float64
        number_of_grids_in_x::Int64
        number_of_grids_in_y::Int64
        Δt::Float64
    end
end

using .Parameters
p = Parameters.Para(
    500.,
    50,
    1.,
    1.,
    80,
    80,
    0.01
    )

function const_diff_mtrx(Nx, Ny, dx, dy)
    """
    2次元空間方向二階差分マトリクスAbigを構築する
    Axは1次元の空間方向二階差分マトリクス
    """
    Abig = zeros(Nx*Ny, Nx*Ny) ;
    Ax = Matrix(-2.0I, Nx, Nx) ;

    for row = 1:Nx
        if row < Nx
            Ax[row, row+1] = 1 ;
        end
        if row > 1
            Ax[row, row-1] = 1 ;
        end
    end
    # boundary condition
    Ax[1, 1] = -1 ;
    Ax[Nx, Nx] = -1 ;

    dd = Matrix(1.0I, Nx, Nx) ;
    for jj=1:Ny
        if jj==1 || jj==Ny 
            @view(Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-1):Nx*jj]) .= Ax ./ dx.^2 .- dd ./ dy.^2 ;
        else
            @view(Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-1):Nx*jj]) .= Ax ./ dx.^2 .- 2 .* dd ./ dy.^2 ;
        end
        if jj!=1
            @view(Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-2):Nx*(jj-1)]) .= dd ./ dy.^2;
        end
        if jj!=Ny
            @view(Abig[1+Nx*(jj-1):Nx*jj, 1+Nx*(jj-0):Nx*(jj+1)]) .= dd ./ dy.^2;
        end
    end
    
    # Abig は特異行列でありこのままでは
    # ポワソン方程式　＋　Neumann境界だと解が一意に定まらないので
    # 1点をu = 0と固定して解とする
    Abig[1,:] .= 0 ;
    Abig[1,1] = 1 ;

    return Abig
end

function main(Para)
    width_of_grid_in_x = Para.length_of_x / Para.number_of_grids_in_x ;   # dx
    width_of_grid_in_y = Para.length_of_y / Para.number_of_grids_in_y ;   # dy

    # define initial physical quantity
    u = zeros(Para.number_of_grids_in_x+1, Para.number_of_grids_in_y+2); # velocity in x direction (u)
    v = zeros(Para.number_of_grids_in_x+2, Para.number_of_grids_in_y+1); # velocity in y direction (v)
    p = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y); # pressure (lagurange multiplier)
    
    # memory allocation
    u_in_cell_center = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);
    u_in_cell_corner = zeros(Para.number_of_grids_in_x+1, Para.number_of_grids_in_y+1);
    v_in_cell_center = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);
    v_in_cell_corner = zeros(Para.number_of_grids_in_x+1, Para.number_of_grids_in_y+1);

    ∇∇ux = zeros(Para.number_of_grids_in_x-1, Para.number_of_grids_in_y);
    ∇∇uy = zeros(Para.number_of_grids_in_x-1, Para.number_of_grids_in_y);
    ∇∇vx = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y-1);
    ∇∇vy = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y-1);

    uuce = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);
    uvco = zeros(Para.number_of_grids_in_x+1, Para.number_of_grids_in_y+1);
    vvce = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);

    Nu = zeros(Para.number_of_grids_in_x-1, Para.number_of_grids_in_y);
    Nv = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y-1);

    b = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);
    RHS = zeros(Para.number_of_grids_in_x * Para.number_of_grids_in_y, 1);
    pp = zeros(Para.number_of_grids_in_x * Para.number_of_grids_in_y, 1);

    velocity_abs = zeros(Para.number_of_grids_in_x, Para.number_of_grids_in_y);

    # construct spatial difference matrix
    A = const_diff_mtrx(Para.number_of_grids_in_x, Para.number_of_grids_in_y, width_of_grid_in_x, width_of_grid_in_y);

    for K = 1:Para.max_time_step
        bctop = 1; # 境界上部の速度 u
        @views(u[:, 1] = - u[:, 2]);    # bottom
        @view(v[:, 1]) .= 0;    # bottom
        @view(u[:,size(u)[2]]) .= 2 .* bctop .- @view(u[:, size(u)[2]-1]);  # top
        @view(v[:, size(v)[2]]) .= 0;   # top
        @view(u[1, :]) .= 0;    # left
        @views(v[1, :] = -v[2, :]); # left
        @view(u[size(u)[1], :]) .= 0;   # right
        @views(v[size(v)[1], :] = -v[size(v)[1]-1, :]);    # right
    
        ∇∇ux .= (@view(u[1:size(u)[1]-2, 2:size(u)[2]-1]) .- 2 .* @view(u[2:size(u)[1]-1, 2:size(u)[2]-1]) .+ @view(u[3:size(u)[1], 2:size(u)[2]-1])) ./ width_of_grid_in_x.^2 ; 
        ∇∇uy .= (@view(u[2:size(u)[1]-1, 1:size(u)[2]-2]) .- 2 .* @view(u[2:size(u)[1]-1, 2:size(u)[2]-1]) .+ @view(u[2:size(u)[1]-1, 3:size(u)[2]])) ./ width_of_grid_in_y.^2 ;
        ∇∇vx .= (@view(v[1:size(v)[1]-2, 2:size(v)[2]-1]) .- 2 .* @view(v[2:size(v)[1]-1, 2:size(v)[2]-1]) .+ @view(v[3:size(v)[1], 2:size(v)[2]-1])) ./ width_of_grid_in_x.^2 ;
        ∇∇vy .= (@view(v[2:size(v)[1]-1, 1:size(v)[2]-2]) .- 2 .* @view(v[2:size(v)[1]-1, 2:size(v)[2]-1]) .+ @view(v[2:size(v)[1]-1, 3:size(v)[2]])) ./ width_of_grid_in_y.^2 ;
    
        # 1. interpolate velocity at cell center/cell cornder
        u_in_cell_center .= (@view(u[1:size(u)[1]-1, 2:size(u)[2]-1]) .+ @view(u[2:size(u)[1], 2:size(u)[2]-1])) ./ 2 ;
        u_in_cell_corner .= (@view(u[:, 1:size(u)[2]-1]) .+ @view(u[:, 2:size(u)[2]])) ./ 2 ;
        v_in_cell_center .= (@view(v[2:size(v)[1]-1, 1:size(v)[2]-1]) .+ @view(v[2:size(v)[1]-1, 2:size(v)[2]])) ./ 2 ;
        v_in_cell_corner .= (@view(v[1:size(v)[1]-1, :]) .+ @view(v[2:size(v)[1], :])) ./ 2 ;
        
        # 2. multiply
        @. uuce = u_in_cell_center * u_in_cell_center ;
        @. uvco = u_in_cell_corner * v_in_cell_corner ;
        @. vvce = v_in_cell_center * v_in_cell_center ;
        
        # 3-1. get derivative for u
        Nu .= (@view(uuce[2:size(uuce)[1], :]) .- @view(uuce[1:size(uuce)[1]-1, :])) ./ width_of_grid_in_x ;
        Nu .+= (@view(uvco[2:size(uvco)[1]-1, 2:size(uvco)[2]]) .- @view(uvco[2:size(uvco)[1]-1, 1:size(uvco)[2]-1])) ./ width_of_grid_in_y ;
    
        # 3-2. get derivative for v
        Nv .= (@view(vvce[:, 2:size(vvce)[2]]) .- @view(vvce[:, 1:size(vvce)[2]-1])) ./ width_of_grid_in_y ;
        Nv .+= (@view(uvco[2:size(uvco)[1], 2:size(uvco)[2]-1]) .- @view(uvco[1:size(uvco)[1]-1, 2:size(uvco)[2]-1])) ./ width_of_grid_in_x ;
    
        # Get intermidiate velocity
        @view(u[2:size(u)[1]-1, 2:size(u)[2]-1]) .= @view(u[2:size(u)[1]-1, 2:size(u)[2]-1]) .+ Para.Δt .* (-Nu .+ (∇∇ux .+ ∇∇uy) ./ Para.Reynolds_number) ;
        @view(v[2:size(v)[1]-1, 2:size(v)[2]-1]) .= @view(v[2:size(v)[1]-1, 2:size(v)[2]-1]) .+ Para.Δt .* (-Nv .+ (∇∇vx .+ ∇∇vy) ./ Para.Reynolds_number) ;
    
        # velocity correction
        # RHS of pressure Poisson eq.
        b .= ((@view(u[2:size(u)[1], 2:size(u)[2]-1]) .- @view(u[1:size(u)[1]-1, 2:size(u)[2]-1])) ./ width_of_grid_in_x .+ (@view(v[2:size(v)[1]-1, 2:size(v)[2]]) .- @view(v[2:size(v)[1]-1, 1:size(v)[2]-1])) ./ width_of_grid_in_y) ;
        RHS .= reshape(b, Para.number_of_grids_in_x * Para.number_of_grids_in_y, 1) ;
        RHS[1] = 0;
        pp .= A\RHS;
        p .= reshape(pp, Para.number_of_grids_in_x, Para.number_of_grids_in_y);
        
        # 圧力勾配を仮の速度場から引いて、新しい速度場
        @view(u[2:size(u)[1]-1, 2:size(u)[2]-1]) .-= (@view(p[2:size(p)[1],:]) .- @view(p[1:size(p)[2]-1,:])) ./ width_of_grid_in_x ;
        @view(v[2:size(v)[1]-1, 2:size(v)[2]-1]) .-= (@view(p[:,2:size(p)[2]]) .- @view(p[:,1:size(p)[2]-1])) ./ width_of_grid_in_y ;
    end

    u_in_cell_center .= (@view(u[1:size(u)[1]-1, 2:size(u)[2]-1]) .+ @view(u[2:size(u)[1], 2:size(u)[2]-1])) ./ 2 ;
    v_in_cell_center .= (@view(v[2:size(v)[1]-1, 1:size(v)[2]-1]) .+ @view(v[2:size(v)[1]-1, 2:size(v)[2]])) ./ 2 ;
        
    velocity_abs .= sqrt.(u_in_cell_center.^2 .+ v_in_cell_center.^2) ;
    contour(velocity_abs, fill = true, c = cgrad(:grays, rev = true)) ;

end

# main(p)
# @time main(p)
@benchmark main(p)