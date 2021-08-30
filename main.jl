# function initgrid(x_g :: Array{Float64, 3}, x_i :: Array{Float64, 2})

# end

function advancetime!(m_i :: Array{Float64, 1}, x_i :: Array{Float64, 2}, v_i :: Array{Float64, 2}, F_i :: Array{Float64, 2}, Δt, box_size)
    n_particles = length(m_i)
    @views for i = 1:n
        r = √(sum(x_i[i,:] .^ 2))
        if r > box_size
            F_i[i,:] = -(1 - box_size/r)*x_i[i,:]
        end
        v_i[i, :] = @. v_i[i, :] + m_i[i]*F_i[i, :]*Δt
        x_i[i, :] = @. x_i[i, :] + v_i[i, :]*Δt
    end
end

function getgriddensity1d(q_i :: Array{Float64, 1}, x_i :: Array{Float64, 2}, N :: Int, box_size)
    grid = 
end
