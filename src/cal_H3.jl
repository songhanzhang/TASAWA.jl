function cal_H3(n_DOF::Int64)
    H3 = zeros(n_DOF,n_DOF)
    H3[[1,n_DOF],[1,n_DOF]] = [1 -1;-1 1]
    for ii = 2:n_DOF-1
        m = ii-1
        H3[ii,ii] = m^2*pi^2/2
    end
    return H3
end
