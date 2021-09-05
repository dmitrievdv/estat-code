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

function initgrid1d(N :: Int, box_size, grid_wide)
    x_g = zeros(N)
    for i = 1:N
        x_g[i] = grid_wide*box_size*((2(i-1)+1)/N - 1)
    end
    return x_g
end

function accumulatechargeanddipole1d!(q_g :: Array{Float64, 1}, D_g :: Array{Float64, 1}, x_g :: Array{Float64, 1}, 
                                        q_i :: Array{Float64, 1}, x_i :: Array{Float64, 1})
    n_particles = length(q_i) 
    N = length(ρ_g)
    Δx_g = x_g[2] - x_g[1]
    grid_size = Δx_g*N/2
    q_g .= 0.0; D_g .= 0.0
    for i = 1:n_particles
        m = Int(floor((grid_size + x_i[i])/Δx_g)) + 1
        if m > N
            m = N
        elseif m < 1
            m = 1
        end
        # println("$m $(x_i[i]) $grid_size $Δx_g")
        x_g = grid_size*((2(m-1)+1)/N - 1)
        q_g[m] += q_i[i]
        D_g[m] += q_i[i]*(x_i[i] - x_g)   
    end
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
        ρ_g[i] += (q_g[1] + (D_g[2] - D_g[1])/Δx_g)*gausscore1d(Δx, σ)
        Δx = x_g[N] - x_g[i]
        ρ_g[i] += (q_g[N] + (D_g[N] - D_g[N-1])/Δx_g)*gausscore1d(Δx, σ)
        for j=2:N-1
            Δx = x_g[j] - x_g[i]
            ρ_g[i] += (q_g[j] + (D_g[j+1] - D_g[j-1])/2Δx_g)*gausscore1d(Δx, σ)
        end
    end
end

n = 10^2
box_size = 1.0
grid_wide = 1.5
x_i = 2*box_size*rand(n) .- box_size
q_i = [ones(n÷2); -ones(n÷2)]

N = 10000
x_g = initgrid1d(N, box_size, grid_wide)
ρ_g = zeros(N); q_g = zeros(N); D_g = zeros(N)

σ = 0.001
accumulatechargeanddipole1d!(q_g, D_g, x_g, q_i, x_i)
griddensity1d!(ρ_g, x_g, q_g, D_g, σ)

ρ_plot = plot(x_g, ρ_g, label = L"\rho_g")

ρ_fft_1 = fft(ρ_g)
ρ_g_1 = real.(ifft(ρ_fft_1))
plot!(ρ_plot, x_g, ρ_g_1, label = "direct IFFT")

# S_g = gausscore1d.(x_g, σ)
q_fft = fft(q_g)
D_fft = fft(D_g)
# S_fft = fft(S_g)
k = fftfreq(N, N/(2*box_size*grid_wide))

# ρ_fft_2 = @. S_fft*(q_fft + D_fft*(1.0im * k * 2π))
# ρ_g_2 = fftshift(real.(ifft(ρ_fft_2)))
# plot!(ρ_plot, x_g, ρ_g_2, label = "IFFT, FFT[S]")



S_fft_th = @. exp(-(2π^2*σ^2*k^2))/(2*grid_wide*box_size/N)
# k = [1:N;]/(2*box_size*grid_wide)
ρ_fft_3 = @. S_fft_th*(q_fft + D_fft*(1.0im * k * 2π))
ρ_g_3 = real.(ifft(ρ_fft_3))
plot!(ρ_plot, x_g, ρ_g_3, label = "IFFT, theorethical FT[S]")


err_plot = plot(k, (@. (abs(ρ_g -ρ_g_1)/n)), label = "direct IFFT" , xlims = (-100, 100))
plot!(err_plot,k, (@. (abs(ρ_g -ρ_g_3)/n)), label = "IFFT, theorethical S(k)")

plot(ρ_plot, err_plot, layout = @layout([A B]), size = (1000, 500), dpi = 100)

ϕ_fft = @. 4π/k^2 * ρ_fft_3
ϕ_fft[1] = 4π/1e-16^2 * ρ_fft_3[1]
ϕ_g = real.(ifft(ϕ_fft))
ϕ_plot = plot(x_g, ϕ_g, label = L"\phi_g")
# vline!(ϕ_plot, x_i[q_i .> 0.0])
# vline!(ϕ_plot, x_i[q_i .< 0.0])

E_fft = @. -1.0im * ρ_fft_3 / (k * 2π) #* S_fft_th
E_fft[1] = -1.0im * ρ_fft_3[1] / (eps(0.0)* 2π) #* S_fft_th[1]
E_g = real.(ifft(E_fft))
E_plot = plot(x_g, E_g, label = L"E_g")
ρ_plot = plot(x_g, ρ_g, label = L"\rho_g")
# vline!(E_plot, x_i[q_i .> 0.0])
# vline!(E_plot, x_i[q_i .< 0.0])

ϕ_plot

plot(ρ_plot, ϕ_plot, E_plot, layout = @layout([A B C]), size = (1500, 500), dpi = 100)