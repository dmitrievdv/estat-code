# function initgrid(x_g :: Array{Float64, 3}, x_i :: Array{Float64, 2})

# end
using Plots
using FFTW
using LaTeXStrings
using Printf

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

function initgrid2d(N :: Int, box_size)
    x_g = zeros(2, N, N)
    Δx_g = 2*box_size/(N-1)
    for i = 1:N, j=1:N
        x_g[1, i, j] = -box_size + (i-1) * Δx_g
        x_g[2, i, j] = -box_size + (j-1) * Δx_g
    end
    return x_g
end

function accumulatechargeanddipole1d!(q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, 
                                        q_i :: Array{Float64, 1}, x_i :: Array{Float64, 1})
    n_particles = length(q_i) 
    N = length(q_g)
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

function accumulatechargeanddipole2d!(q_g :: Array{Float64, 2}, D_g :: Array{Float64, 3}, x_g :: Array{Float64, 3}, 
                                        q_i :: Array{Float64, 1}, x_i :: Array{Float64, 2})
    n_particles = length(q_i) 
    N = size(q_g)[1]
    Δx_g = x_g[1,2,1] - x_g[1,1,1]
    grid_size = Δx_g*(N-1)/2
    q_g .= 0.0; D_g .= 0.0
    for i = 1:n_particles
        m_x = Int(floor(2*(x_i[1,i]+1.0)/Δx_g))
        m_x = ((m_x+1)÷2) + 1

        m_y = Int(floor(2*(x_i[2,i]+1.0)/Δx_g))
        m_y = ((m_y+1)÷2) + 1

        x_g = -grid_size + (m_x-1)*Δx_g
        y_g = -grid_size + (m_y-1)*Δx_g
        q_g[m_x, m_y] += q_i[i]
        D_g[m_x, m_y, 1] += q_i[i]*(x_i[1, i] - x_g)
        D_g[m_x, m_y, 2] += q_i[i]*(x_i[2, i] - y_g)      
    end
    q_g[1,:] = q_g[1,:] .+ q_g[N,:]; q_g[N,:] .= q_g[1,:]
    q_g[:,1] = q_g[:,1] .+ q_g[:,N]; q_g[:,N] .= q_g[:,1]

    D_g[1,:,:] = D_g[1,:,:] .+ D_g[N,:,:]; D_g[N,:,:] .= D_g[1,:,:]
    D_g[:,1,:] = D_g[:,1,:] .+ D_g[:,N,:]; D_g[:,N,:] .= D_g[:,1,:]
end

function gausscore1d(Δx, σ)
    return exp(-0.5*(Δx/σ)^2)/√(2π)/σ
end

function gausscore2d(Δx :: Vector{Float64}, σ)
    Δx² = sum(Δx .^ 2)
    return exp(-0.5*Δx²/(σ)^2)/2π/σ^2
end

function griddensity1d!(ρ_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, σ :: Real)
    N = length(x_g)
    Δx_g = x_g[2] - x_g[1]
    for i = 1:N
        ρ_g[i] = 0.0
        Δx = x_g[1] - x_g[i]
        ρ_g[i] += (q_g[1] + (D_g[N-1] - D_g[2])/2Δx_g)*gausscore1d(Δx, σ)
        Δx = x_g[N] - x_g[i]
        ρ_g[i] += (q_g[N] + (D_g[2] - D_g[N-1])/2Δx_g)*gausscore1d(Δx, σ)
        for j=2:N-1
            Δx = x_g[j] - x_g[i]
            ρ_g[i] += (q_g[j] + (D_g[j+1] - D_g[j-1])/2Δx_g)*gausscore1d(Δx, σ)
        end
    end
end

function griddensity2d!(ρ_g :: Array{Float64, 2}, x_g :: Array{Float64, 3}, q_g :: Array{Float64, 2}, D_g :: Array{Float64, 3}, σ :: Real)
    N = size(ρ_g)[1]
    Δx_g = x_g[1,2,1] - x_g[1,1,1]
    # println(Δx_g)
    Δx = zeros(2)
    @views for i = 1:N, j=1:N
        ρ_g[i,j] = 0.0

        # --- corner points -------------
        Δx .= x_g[:,1,1] .- x_g[:,i,j]
        ρ_g[i,j] += (q_g[1,1] + (D_g[N-1,1,1] - D_g[2,1,1] + D_g[1,N-1,2] - D_g[1,2,2])/2Δx_g)*gausscore2d(Δx, σ)
        Δx .= x_g[:,1,N] .- x_g[:,i,j]
        ρ_g[i,j] += (q_g[1,N] + (D_g[N-1,N,1] - D_g[2,N,1] + D_g[1,N-1,2] - D_g[1,2,2])/2Δx_g)*gausscore2d(Δx, σ)
        Δx .= x_g[:,N,1] .- x_g[:,i,j]
        ρ_g[i,j] += (q_g[N,1] + (D_g[N-1,1,1] - D_g[2,1,1] + D_g[N,N-1,2] - D_g[N,2,2])/2Δx_g)*gausscore2d(Δx, σ)
        Δx .= x_g[:,N,N] .- x_g[:,i,j]
        ρ_g[i,j] += (q_g[N,N] + (D_g[N-1,N,1] - D_g[2,N,1] + D_g[N,N-1,2] - D_g[N,2,2])/2Δx_g)*gausscore2d(Δx, σ)
        # -------------------------------

        # --- left and right border -----
        for j‵ = 2:N-1
            Δx .= x_g[:,1,j‵] .- x_g[:,i,j]
            ρ_g[i,j] += (q_g[1,j‵] + (D_g[N-1,j‵,1] - D_g[2,j‵,1] + D_g[1,j‵+1,2] - D_g[1,j‵-1,2])/2Δx_g)*gausscore2d(Δx, σ)

            Δx .= x_g[:,N,j‵] .- x_g[:,i,j]
            ρ_g[i,j] += (q_g[N,j‵] + (D_g[N-1,j‵,1] - D_g[2,j‵,1] + D_g[N,j‵+1,2] - D_g[N,j‵-1,2])/2Δx_g)*gausscore2d(Δx, σ)
        end
        # --------------------------------

        # --- top and bottom border ------
        for i‵ = 2:N-1
            Δx .= x_g[:,i‵,1] .- x_g[:,i,j]
            ρ_g[i,j] += (q_g[i‵,1] + (D_g[i‵,N-1,2] - D_g[i‵,2,2] + D_g[i‵+1,1,1] - D_g[i‵-1,1,1])/2Δx_g)*gausscore2d(Δx, σ)

            Δx .= x_g[:,i‵,N] .- x_g[:,i,j]
            ρ_g[i,j] += (q_g[i‵,N] + (D_g[i‵,N-1,2] - D_g[i‵,2,2] + D_g[i‵+1,N,1] - D_g[i‵-1,N,1])/2Δx_g)*gausscore2d(Δx, σ)
        end
        # --------------------------------

        # --- everything else ------------
        for i‵=2:N-1, j‵=2:N-1
            Δx .= x_g[:,i‵,j‵] .- x_g[:,i,j]
            ρ_g[i,j] += (q_g[i‵,j‵] + (D_g[i‵,j‵+1,2] - D_g[i‵,j‵-1,2] + D_g[i‵+1,N,1] - D_g[i‵-1,N,1])/2Δx_g)*gausscore2d(Δx, σ)
        end
        # --------------------------------
    end
end

"""
    poissonequation1d(n :: Int, N :: Int, box_size :: Real, σ :: Real)

Solve 1d poisson equation for `n` positively charged and `n` negatevely charged `σ`-sized particles
randomly distributed in `[-box_size:box_size]` using FFT on grid with `N` points.
"""
function poissonequation1d(n_p :: Int, N :: Int, box_size, σ)
    n = 2n_p # must be even
 
    x_i = 2*box_size*rand(n) .- box_size # n random particles in the box
    q_i = [ones(n÷2); -ones(n÷2)] # half have charge +1 and half -- -1
    # q_i = [rand(-1.0:2:1) for i =1:n] # --- for some reason this doesn't work properly

    x_g = initgrid1d(N, box_size) # grid creation
    ρ_g = zeros(N); q_g = zeros(N); D_g = zeros(N) # grid values: density, charge and dipole moment

    ρ_plot = plot()

    accumulatechargeanddipole1d!(q_g, D_g, x_g, q_i, x_i) # accumulation of charge and dipole moment on the grid

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
end

function poissonequation2d(n_p :: Int, N :: Int, box_size, σ)
    n = 2n_p # must be even
 
    x_i = 2*box_size*rand(2,n) .- box_size # n random particles in the box
    q_i = [ones(n÷2); -ones(n÷2)] # half have charge +1 and half -- -1
    # q_i = [rand(-1.0:2:1) for i =1:n] # --- for some reason this doesn't work properly

    x_g = initgrid2d(N, box_size) # grid creation
    ρ_g = zeros(N,N); q_g = zeros(N,N); D_g = zeros(N,N,2) # grid values: density, charge and dipole moment

    # ρ_plot = plot()

    accumulatechargeanddipole2d!(q_g, D_g, x_g, q_i, x_i) # accumulation of charge and dipole moment on the grid
    
    
    # griddensity2d!(ρ_g, x_g, q_g, D_g, σ) # calculation of density on grid (don't need this if you trust in FFT)
    # #                                       # (but it wouldn't properly work near the edges)
    
    # ρ_g
    # # ρ_plot = plot(x_g, ρ_g, label = L"\rho_g") # uncomment this and the line above (griddensity1d...) if you don't trust in FFT

    # ----- calculation of ρ_g using FFT ---------
    q_fft = fft(q_g)
    D_x_fft = fft(D_g[:,:,1])
    D_y_fft = fft(D_g[:,:,2])
    k_x = fftfreq(N, N/(2*box_size))
    k_y = fftfreq(N, N/(2*box_size))
    # S_fft_th = @. exp(-(2π^2*σ^2*k^2))/(2*box_size/N)^2
    ρ_fft_3 = zeros(ComplexF64,N,N)
    E_fft = zeros(ComplexF64,N,N)
    ϕ_fft = zeros(ComplexF64,N,N)
    for i = 1:N, j=1:N
        S_k = exp(-(2π^2*σ^2*(k_x[i]^2 + k_y[j]^2)))/(2*box_size/N)^2
        ρ_fft_3[i,j] = S_k*(q_fft[i,j] + D_x_fft[i,j]*(1.0im * k_x[i] * 2π) + D_y_fft[i,j]*(1.0im * k_y[j] * 2π))
        ϕ_fft[i,j] = @. 4π/((k_x[i]^2 + k_y[j]^2) * (2π)^2) * ρ_fft_3[i,j] 
    end
    ϕ_fft[1,1] = 0.0 + 0.0im
    # ρ_fft_3 = @. S_fft_th*(q_fft + D_fft*(1.0im * k * 2π))
    ρ_g_3 = real.(ifft(ρ_fft_3))
    # ---------------------------------------------

    # --- solving poisson equation ---------------------------------
    ϕ_g = real.(ifft(ϕ_fft))
    # ϕ_plot = plot(x_g, ϕ_g, label = L"\phi_g")
    # --------------------------------------------------------------

    
    out = open("2dgrid.dat", "w")

    for i=1:N
        for j=1:N
            @printf(out, "%9.5f %9.5f %15.5e %15.5e\n", x_g[1,i,j], x_g[2,i,j], ρ_g_3[i,j], ϕ_g[i,j])
        end
        println(out, "")
    end

    close(out)

    out = open("2dparticles.dat", "w")

    for i=1:n
        @printf(out, "%9.5f %9.5f %9.5f\n", x_i[1,i], x_i[2,i], q_i[i])
    end

    close(out)

    # xs = x_g[1,:,1]
    # ys = x_g[2,1,:]
    # ρ_plot = heatmap(ys, xs, ρ_g_3, aspect_ratio = :equal)
    # ϕ_plot = heatmap(ys, xs, ϕ_g, aspect_ratio = :equal)
    # positives = q_i .> 0.0
    # negatives = q_i .< 0.0
    # scatter!(ρ_plot, x_i[2,positives], x_i[1,positives], mc = :red)
    # scatter!(ρ_plot, x_i[2,negatives], x_i[1,negatives], mc = :blue)
    # scatter!(ϕ_plot, x_i[2,positives], x_i[1,positives], mc = :red)
    # scatter!(ϕ_plot, x_i[2,negatives], x_i[1,negatives], mc = :blue)

    # return ρ_g, ϕ_g

    # plot(ρ_plot, ϕ_plot, layout = @layout([A B]), size = (1000, 500), dpi = 100)

    # # --- finding E field or force on particle ---------------------
    # E_fft = @. -1.0im * ρ_fft_3 / (k * 2π) #* S_fft_th 
    # E_fft[1] = 0.0# -1.0im * ρ_fft_3[1] / (eps(0.0)* 2π) #* S_fft_th[1] --- don't need this k
    # E_g = real.(ifft(E_fft))
    # E_plot = plot(x_g, E_g, label = L"E_g")
    # # --------------------------------------------------------------


    # # ------ uncomment vlines if you want to add vertical lines on charges
    # # vline!(ρ_plot, x_i[q_i .> 0.0], lc = :red)
    # # vline!(ρ_plot, x_i[q_i .< 0.0], lc = :blue)

    # # vline!(ϕ_plot, x_i[q_i .> 0.0], lc = :red)
    # # vline!(ϕ_plot, x_i[q_i .< 0.0], lc = :blue)

    # # vline!(E_plot, x_i[q_i .> 0.0], lc = :red)
    # # vline!(E_plot, x_i[q_i .< 0.0], lc = :blue)
    # # ---------------------------------------------------------------------


    # plot(ρ_plot, ϕ_plot, E_plot, layout = @layout([A B C]), size = (1500, 500), dpi = 100, xlims = (-1.0, 1.0))
end

@time poissonequation2d(10^5,500,1.0,0.01)

# @time poissonequation1d(1000, 1000, 1.0, 0.01)