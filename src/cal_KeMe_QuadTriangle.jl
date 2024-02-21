function cal_KeMe_QuadraticTriangular(x,y,E,ρ,ν)

    #= Benchmark Example
    x1 = 0
    y1 = 0
    x2 = 1.5
    y2 = 0.2
    x3 = 0.6
    y3 = 0.9

    x4 = (x1+x2)/2
    y4 = (y1+y2)/2
    x5 = (x2+x3)/2
    y5 = (y2+y3)/2
    x6 = (x1+x3)/2
    y6 = (y1+y3)/2

    x = [ x1 x2 x3 x4 x5 x6 ]'
    y = [ y1 y2 y3 y4 y5 y6 ]'

    E = 2e11
    ρ  = 7850
    ν = 0.3
    =#

    λ = E*ν/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))

    c11 = λ+2*μ
    c22 = λ+2*μ
    c33 = λ+2*μ
    c12 = λ
    c13 = λ
    c23 = λ
    c44 = μ

    C = [ c11 c12 c13 0   0   0
          c12 c22 c23 0   0   0
          c13 c23 c33 0   0   0
          0   0   0   c44 0   0
          0   0   0   0   c44 0
          0   0   0   0   0   c44 ]

    L = [ 1 0 0 0 0 0 0 0 0
          0 0 0 0 1 0 0 0 0
          0 0 0 0 0 0 0 0 1
          0 1 0 1 0 0 0 0 0
          0 0 1 0 0 0 1 0 0
          0 0 0 0 0 1 0 1 0 ]

    K1e = zeros(18,18)
    K2e = zeros(18,18)
    K3e = zeros(18,18)
    Me  = zeros(18,18)

    Gauss = [ 0.0915762135  0.8168475730  0.1099517437
              0.0915762135  0.0915762135  0.1099517437
              0.8168475730  0.0915762135  0.1099517437
              0.4459484909  0.1081030182  0.2233815897
              0.4459484909  0.4459484909  0.2233815897
              0.1081030182  0.4459484909  0.2233815897 ]

    for i_Gauss = 1:6
        ξ = Gauss[i_Gauss,1]
        η = Gauss[i_Gauss,2]
        H = Gauss[i_Gauss,3]
        Nb = zeros(6)
        Nb[1] = (1-ξ-η) * (1-2*ξ-2*η)
        Nb[2] = ξ * (2*ξ-1)
        Nb[3] = η * (2*η-1)
        Nb[4] = 4 * ξ * (1-ξ-η)
        Nb[5] = 4 * ξ * η
        Nb[6] = 4 * η * (1-ξ-η)
        dN_dξ = zeros(6)
        dN_dξ[1] = 4*η + 4*ξ - 3
        dN_dξ[2] = 4*ξ - 1
        dN_dξ[3] = 0
        dN_dξ[4] = 4 - 8*ξ - 4*η
        dN_dξ[5] = 4*η
        dN_dξ[6] = -4*η
        dN_dη = zeros(6)
        dN_dη[1] = 4*η + 4*ξ - 3
        dN_dη[2] = 0
        dN_dη[3] = 4*η - 1
        dN_dη[4] = -4*ξ
        dN_dη[5] = 4*ξ
        dN_dη[6] = 4 - 4*ξ - 8*η
        dx_dξ = dN_dξ' * x
        dy_dξ = dN_dξ' * y
        dx_dη = dN_dη' * x
        dy_dη = dN_dη' * y
        J = [ dx_dξ  dy_dξ
              dx_dη  dy_dη ]
        dN_dx = zeros(6)
        dN_dy = zeros(6)
        for ii = 1:6
            Temp = J \ [ dN_dξ[ii]; dN_dη[ii] ]
            dN_dx[ii] = Temp[1,1]
            dN_dy[ii] = Temp[2,1]
        end
        B1 = [ zeros(6,18)
               kron(Nb',    [ 1 0 0; 0 1 0; 0 0 1 ]) ]
        B2 = [ kron(dN_dx', [ 1 0 0; 0 1 0; 0 0 1 ])
               kron(dN_dy', [ 1 0 0; 0 1 0; 0 0 1 ])
               zeros(3,18)]
        N = kron(Nb', [ 1 0 0; 0 1 0; 0 0 1 ])
        K1e += B1'*L'*C*L*B1*abs(det(J))*H
        K2e += (B2'*L'*C*L*B1 - B1'*L'*C*L*B2)*abs(det(J))*H
        K3e += B2'*L'*C*L*B2*abs(det(J))*H
        Me += ρ*transpose(N)*N*abs(det(J))*H
        if det(J) < 0
            println(string("\nAttension: Negative Jacobian for element ", i_e))
        end
    end

    return K1e, K2e, K3e, Me

end
