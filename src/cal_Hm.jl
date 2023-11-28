function cal_Hm(n_DOF::Int64)
    Hm = zeros(n_DOF,n_DOF)
    # Block (1,1)
    ii = 1
    jj = 1
    Hm[ii,jj] = 1/3
    # Block (1,2)
    ii = 1
    for jj = 2:n_DOF-1
        n = jj-1
        Hm[ii,jj] = 1/(n*pi)
    end
    # Block (1,3)
    ii = 1
    jj = n_DOF
    Hm[ii,jj] = 1/6
    # Block (2,1)
    jj = 1
    for ii = 2:n_DOF-1
        m = ii-1
        Hm[ii,jj] = 1/(m*pi)
    end
    # Block(2,2)
    for ii = 2:n_DOF-1
        for jj = 2:n_DOF-1
            m = ii-1
            n = jj-1
            Hm[ii,jj] = 1/2*(m == n)
        end
    end
    # Block(2,3)
    jj = n_DOF
    for ii = 2:n_DOF-1
        m = ii-1
        Hm[ii,jj] = (-1)^(m+1)/(m*pi)
    end
    # Block(3,1)
    ii = n_DOF
    jj = 1
    Hm[ii,jj] = 1/6
    # Block(3,2)
    ii = n_DOF
    for jj = 2:n_DOF-1
        n = jj-1
        Hm[ii,jj] = (-1)^(n+1)/(n*pi)
    end
    # Block(3,3)
    ii = n_DOF;
    jj = n_DOF;
    Hm[ii,jj] = 1/3
    return Hm
end
