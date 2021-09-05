# function initgrid(x_g :: Array{Float64, 3}, x_i :: Array{Float64, 2})

# end
using Plots
using FFTW
using LaTeXStrings

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
    Δx_g = 2*box_size/(N-1)
    for i = 1:N
        x_g[i] = -box_size + (i-1) * Δx_g
    end
    return x_g
end

function accumulatechargeanddipole1d!(q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, 
                                        q_i :: Array{Float64, 1}, x_i :: Array{Float64, 1})
    n_particles = length(q_i) 
    N = length(ρ_g)
    Δx_g = x_g[2] - x_g[1]
    grid_size = Δx_g*(N-1)/2
    q_g .= 0.0; D_g .= 0.0
    for i = 1:n_particles
        m = Int(floor(2*(x_i[i]+1.0)/Δx_g))
        # print("$(x_i[i]) $m ")
        m = ((m+1)÷2) + 1
        # println("$m")
        # println("$m $(x_i[i]) $grid_size $Δx_g")
        x_g = -grid_size + (m-1)*Δx_g
        q_g[m] += q_i[i]
        D_g[m] += q_i[i]*(x_i[i] - x_g)   
    end
    q_g[1] += q_g[N]; q_g[N] = q_g[1]
    D_g[1] += D_g[N]; D_g[N] = D_g[1]
end

function gausscore1d(Δx, σ)
    return exp(-0.5*(Δx/σ)^2)/√(2π)/σ
end

function griddensity1d!(ρ_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, σ :: Real)
    N = length(x_g)
    Δx_g = x_g[2] - x_g[1]
    for i = 1:N
        ρ_g[i] = 0.0
        Δx = x_g[1] - x_g[i]
        ρ_g[i] += (q_g[1] + (D_g[N-1] - D_g[1])/2Δx_g)*gausscore1d(Δx, σ)
        Δx = x_g[N] - x_g[i]
        ρ_g[i] += (q_g[N] + (D_g[2] - D_g[N-1])/2Δx_g)*gausscore1d(Δx, σ)
        for j=2:N-1
            Δx = x_g[j] - x_g[i]
            ρ_g[i] += (q_g[j] + (D_g[j+1] - D_g[j-1])/2Δx_g)*gausscore1d(Δx, σ)
        end
    end
end

n = 10^2 # number of particles
n = 2n # must be even
box_size = 1.0 # grid from -box_size to +box_size
 
x_i = 2*box_size*rand(n) .- box_size # n random particles in the box
q_i = [ones(n÷2); -ones(n÷2)] # half have charge +1 and half -- -1
# q_i = [rand(-1.0:2:1) for i =1:n] # --- for some reason this doesn't work properly

N = 1000 # number of grid points
x_g = initgrid1d(N, box_size) # grid creation
ρ_g = zeros(N); q_g = zeros(N); D_g = zeros(N) # grid values: density, charge and dipole moment

ρ_plot = plot()

accumulatechargeanddipole1d!(q_g, D_g, x_g, q_i, x_i) # accumulation of charge and dipole moment on the grid

σ = 0.01 # gaussian width --- each particle is a gaussian (see gausscore1d)
# griddensity1d!(ρ_g, x_g, q_g, D_g, σ) # calculation of density on grid (don't need this if you trust in FFT)
#                                       # (but it wouldn't properly work near the edges)
# ρ_plot = plot(x_g, ρ_g, label = L"\rho_g") # uncomment this and the line above (griddensity1d...) if you don't trust in FFT

# ----- calculation of ρ_g using FFT ---------
q_fft = fft(q_g)
D_fft = fft(D_g)
k = fftfreq(N, N/(2*box_size))
S_fft_th = @. exp(-(2π^2*σ^2*k^2))/(2*box_size/N)
ρ_fft_3 = @. S_fft_th*(q_fft + D_fft*(1.0im * k * 2π))
ρ_g_3 = real.(ifft(ρ_fft_3))
plot!(ρ_plot, x_g, ρ_g_3, label = L"\rho_g,\ \mathrm{FFT}")
# ---------------------------------------------

# --- solving poisson equation ---------------------------------
ϕ_fft = @. 4π/(k * 2π)^2 * ρ_fft_3 
ϕ_fft[1] = 0.0#  4π/(1e-20) * ρ_fft_3[1] --- don't need this k
ϕ_g = real.(ifft(ϕ_fft))
ϕ_plot = plot(x_g, ϕ_g, label = L"\phi_g")
# --------------------------------------------------------------

# --- finding E field or force on particle ---------------------
E_fft = @. -1.0im * ρ_fft_3 / (k * 2π) #* S_fft_th 
E_fft[1] = 0.0# -1.0im * ρ_fft_3[1] / (eps(0.0)* 2π) #* S_fft_th[1] --- don't need this k
E_g = real.(ifft(E_fft))
E_plot = plot(x_g, E_g, label = L"E_g")
# --------------------------------------------------------------


# ------ uncomment vlines if you want to add vertical lines on charges
# vline!(ρ_plot, x_i[q_i .> 0.0], lc = :red)
# vline!(ρ_plot, x_i[q_i .< 0.0], lc = :blue)

# vline!(ϕ_plot, x_i[q_i .> 0.0], lc = :red)
# vline!(ϕ_plot, x_i[q_i .< 0.0], lc = :blue)

# vline!(E_plot, x_i[q_i .> 0.0], lc = :red)
# vline!(E_plot, x_i[q_i .< 0.0], lc = :blue)
# ---------------------------------------------------------------------


plot(ρ_plot, ϕ_plot, E_plot, layout = @layout([A B C]), size = (1500, 500), dpi = 100, xlims = (-1.0, 1.0))