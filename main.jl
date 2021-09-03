# function initgrid(x_g :: Array{Float64, 3}, x_i :: Array{Float64, 2})

# end

function advancetime!(m_i :: Array{Float64, 1}, x_i :: Array{Float64, 2}, v_i :: Array{Float64, 2}, F_i :: Array{Float64, 2}, Δt, box_size; stiffness = 10.0)
    n_particles = length(m_i)
    @views for i = 1:n_particles
        r = √(sum(x_i[i,:] .^ 2))
        if r > box_size
            F_i[i,:] = -stiffness*(1 - box_size/r)*x_i[i,:]
        end
        v_i[i, :] = @. v_i[i, :] + m_i[i]*F_i[i, :]*Δt
        x_i[i, :] = @. x_i[i, :] + v_i[i, :]*Δt
    end
end

function initgrid1d(N :: Int, box_size)
    x_g = zeros(N)
    for i = 1:N
        x_g[i] = box_size*((2(i-1)+1)/N - 1)
    end
end

function accumulatechargeanddipole1d!(q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, 
                                        q_i :: Array{Float64, 1}, x_i :: Array{Float64, 1}, N :: Int, box_size)
    n_particles = length(q_i) 
    N = length(ρ_g)
    Δx = 2*box_size/N
    q_g .= 0.0; D_g .= 0.0
    for i = 1:n_particles
        m = Int(floor((box_size + x_i[i])/Δx)) + 1
        if m > N
            m = N
        elseif m < 1
            m = 1
        end
        x_g = box_size*((2(m-1)+1)/N - 1)
        q_g[m] += q_i[i]
        D_g[m] += x_i[i] - x_g[m]   
    end
end

function griddensity1d!(ρ_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, σ :: Real)
    N = length(x_g)
    for i = 1:N
        σ/√(2π)*exp(-(x))
    end
end
