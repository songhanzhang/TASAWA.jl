function cal_type(var,cr)
    if abs(real(var))>cr && abs(imag(var))<cr
        var_type = "real"
    elseif abs(real(var))<cr && abs(imag(var))>cr
        var_type = "imag"
    elseif abs(real(var))>cr && abs(imag(var))>cr
        var_type = "cplx"
    else abs(real(var))<cr && abs(real(var))<cr
        var_type = "zero"
    end
    return var_type
end
