function div_k_types(k_save; tol = 1e-6)

    k_real = zeros(size(k_save))*NaN
    k_imag = zeros(size(k_save))*NaN
    k_cplx = zeros(size(k_save))*(NaN+NaN*1im)
    for ii = 1:size(k_save,1)
        for jj = 1:size(k_save,2)
            k = k_save[ii,jj]
            if abs(real(k))>1e-3 && abs(imag(k))<tol
                k_real[ii,jj] = real(k_save[ii,jj])
                # cg_real[ii,jj] = real(cg_save[ii,jj])
            elseif abs(real(k))<1e-3 && abs(imag(k))>tol
                k_imag[ii,jj] = imag(k_save[ii,jj])
            elseif abs(real(k))>1e-3 && abs(imag(k))>tol
                k_cplx[ii,jj] = k_save[ii,jj]
            end
        end
    end

    return k_real, k_imag, k_cplx

end
