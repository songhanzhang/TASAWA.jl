function cal_H2(n_DOF::Int64)
    H2 = zeros(n_DOF,n_DOF)
    # Block (1,1)
    ii = 1
    jj = 1
    H2[ii,jj] = -1/2
    # Block (1,2)
    ii = 1
    for jj = 2:n_DOF-1
        n = jj-1
        H2[ii,jj] = 2/(n*pi)*(mod(n,2) == 1)
    end
    # Block (1,3)
    ii = 1
    jj = n_DOF
    H2[ii,jj] = 1/2
    # Block (2,1)
    jj = 1
    for ii = 2:n_DOF-1
        m = ii-1
        H2[ii,jj] = -2/(m*pi)*(mod(m,2) == 1)
    end
    # Block(2,2)
    for ii = 2:n_DOF-1
        for jj = 2:n_DOF-1
            m = ii-1
            n = jj-1
            if m != n
                H2[ii,jj] = -(m*n*((-1)^(m + n) - 1))/(m^2 - n^2)
            end
        end
    end
    # Block(2,3)
    jj = n_DOF
    for ii = 2:n_DOF-1
        m = ii-1
        H2[ii,jj] = (2*sin((pi*m)/2)^2)/(m*pi)*(mod(m,2) == 1)
    end
    # Block(3,1)
    ii = n_DOF
    jj = 1
    H2[ii,jj] = -1/2
    # Block(3,2)
    ii = n_DOF
    for jj = 2:n_DOF-1
        n = jj-1
        H2[ii,jj] = -2/(n*pi)*(mod(n,2) == 1)
    end
    # Block(3,3)
    ii = n_DOF
    jj = n_DOF
    H2[ii,jj] = 1/2
    return H2
end
