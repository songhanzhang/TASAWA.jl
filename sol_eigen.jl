function sol_eigen(K1g,K2g_hat,K3g,f_ax)

    n_DOF = size(K1g,1)

    k_save = zeros(2*n_DOF,length(f_ax))*(0+0im)
    Ψ_save = zeros(2*n_DOF, 2*n_DOF, length(f_ax))*(0+0im)

    for (i_freq,freq) in enumerate(f_ax)
        @printf "Current progress: %.2f%%\n" i_freq/length(f_ax)*100
        ω = 2*pi*freq
        A_eig = [zeros(n_DOF,n_DOF)  K3g-ω^2*Mg
                K3g-ω^2*Mg        K2g_hat]
        B_eig = [K3g-ω^2*Mg        zeros(n_DOF,n_DOF)
                zeros(n_DOF,n_DOF)  -K1g]
        (Λ,Φ) = eigen(A_eig,B_eig)
        k_save[:,i_freq] = Λ
        Ψ_save[:,:,i_freq] = Φ
    end

    return k_save, Ψ_save

end