using MAT
using LinearAlgebra
using SparseArrays
using Plots
using Printf
using Measures

include("/Users/songhan.zhang/Documents/GitHub/TASAWA.jl/src/TASAWA.jl")

#=
file = matopen("/Users/songhan.zhang/Documents/MATLAB/2023-QuadraticTriangle/mesh_circle.mat")
Nodes_mesh = read(file, "Nodes_mesh")
Elements_mesh = read(file, "Elements_mesh")
n_nodes = size(Nodes_mesh,1)
n_elements = size(Elements_mesh,1)

Nodes = [ 1:1:n_nodes  Nodes_mesh ]
Elements = Array{Any,2}(undef,n_elements,5)
for i_e = 1:n_elements
    Elements[i_e,1] = i_e
    Elements[i_e,2] = "2D_QuadTriangle"
    Elements[i_e,3] = 1
    Elements[i_e,4] = 1
    Elements[i_e,5] = Elements_mesh[i_e,:]
end
=#

work_path = "/Users/songhan.zhang/Documents/Julia/2023-Julia-v1205-SAFE/"

Nodes = [
    1  0.00  0.00
    2  0.05  0.00
    3  0.05*cos(pi/4)  0.05*sin(pi/4)
    4  0.05*cos(2*pi/4)  0.05*sin(2*pi/4)
    5  0.05*cos(3*pi/4)  0.05*sin(3*pi/4)
    6  0.05*cos(4*pi/4)  0.05*sin(4*pi/4)
    7  0.05*cos(5*pi/4)  0.05*sin(5*pi/4)
    8  0.05*cos(6*pi/4)  0.05*sin(6*pi/4)
    9  0.05*cos(7*pi/4)  0.05*sin(7*pi/4)
    10   0.025   0.000
    11   0.000   0.025
    12  -0.025   0.000
    13   0.000  -0.025
]
n_nodes = size(Nodes,1)

Elements = [
    1  "2D_QuadTriangle"  1  1  (1, 2, 4, 10, 3, 11)
    # 2  "2D_QuadTriangle"  1  1  (1, 4, 6, 11, 5, 12)
    # 3  "2D_QuadTriangle"  1  1  (1, 6, 8, 12, 7, 13)
    # 4  "2D_QuadTriangle"  1  1  (1, 8, 2, 13, 9, 10)
]
n_elements = size(Elements,1)

Materials = [ 1  (2e11, 7800, 0.33) ]
Reals = [ 1  (1) ]

fig_model = plot(
    size = (600,600),
    dpi = 900,
    legend = false,
    grid = false,
    frame_style = :box,
    tickfontsize = 10,
    aspect_ratio = :equal
)
plot_model_elements(Nodes, Elements)
plot_model_nodes(Nodes)
savefig(string(work_path, "fig_model.pdf"))

(n_DOF, list_DOF) = cal_list_DOF(n_nodes, [1,2,3])

(K1g, K2g, K3g, Mg) = cal_KgMg(Nodes, Elements, Materials, list_DOF)

T = Matrix(I,n_DOF,n_DOF)*(1+0im)
for i_DOF = 1:floor(Int,n_DOF/3)
    T[(i_DOF-1)*3+3,(i_DOF-1)*3+3] = 1im
end
K2g_hat = real(1im*T'*K2g*T)
f_ax = (0.00001:0.00001:0.06)*1e6
(k_save,Ψ_save) = sol_eigen(K1g,K2g_hat,K3g,f_ax)
(k_real, k_imag, k_cplx) = div_k_types(k_save, tol = 1e-3)

fig_kw = plot(
    size = (600,700),
    dpi = 900,
    grid = false,
    legend = false,
    frame_style = :box,
    tickfontsize = 10
)
for ii = 1:n_DOF*2
    scatter!(
        -abs.(k_imag[ii,:]), 2*pi*f_ax/1e6, label = "",
        markersize = 0.5, markerstrokecolor = :green, color = :green
    )
    scatter!(
        abs.(real(k_cplx[ii,:])), 2*pi*f_ax/1e6, label = "",
        markersize = 0.5, markerstrokecolor = :dodgerblue, color = :dodgerblue
    )
    scatter!(
        abs.(-imag(k_cplx[ii,:])), 2*pi*f_ax/1e6, label = "",
        markersize = 0.5, markerstrokecolor = :dodgerblue, color = :dodgerblue
    )
    scatter!(
        abs.(k_real[ii,:]), 2*pi*f_ax/1e6, label = "",
        markersize = 1.0, markerstrokecolor = :hotpink1, color = :hotpink1
    )
end
plot!(
    [0,0], [0,0.40], label = "",
    w = 1.5, linecolor = :black, linestyle = :dash
)
plot!(xlims = (-30,60), ylims = (0,0.35))
xlabel!("k (m⁻¹)")
ylabel!("ω (10⁶ rad/s)")
savefig(string(work_path, "fig_kw.pdf"))

# plot!(xlims = (22.5,25.5), ylims = (0.272,0.273))

i_freq = 4340
i_mode = 60
ω_sel = f_ax[i_freq]*2*pi
k_save[i_mode,i_f]

Q = Ψ_save[1:n_DOF,i_mode,i_freq]

m = real(transpose(Q)*Mg*Q)
k1 = real(transpose(Q)*K1g*Q)
k2 = real(transpose(Q)*K2g_hat*Q)
k3 = real(transpose(Q)*K3g*Q)

k_peak = -k2/2/k1

k = k_save[i_mode,i_freq]
ω_peak = sqrt((k_peak^2*k1 + k_peak*k2 + k3)/(m))/1e6

ω = f_ax[i_freq]*2*pi
(-ω^2*Mg + K3g + k*K2g_hat + k^2*K1g)*Q
