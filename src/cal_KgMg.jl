function cal_KgMg(Nodes, Elements, Materials, list_DOF)

    n_gDOF = size(list_DOF,1)
    n_elements = size(Elements,1)

    K1g = zeros(n_gDOF,n_gDOF)
    K2g = zeros(n_gDOF,n_gDOF)
    K3g = zeros(n_gDOF,n_gDOF)
    Mg  = zeros(n_gDOF,n_gDOF)

    for i_e = 1:n_elements
        if Elements[i_e,2] == "1D_GDSA"
            node_1 = Elements[i_e,5][1]
            node_2 = Elements[i_e,5][2]
            n_mode = Elements[i_e,5][3]
            n_bDOF = 2 + n_mode
            Hm = cal_Hm(n_bDOF)
            H2 = cal_H2(n_bDOF)
            H3 = cal_H3(n_bDOF)
            i_mat = Elements[i_e,3]
            rho = Materials[i_mat,2]
            E = Materials[i_mat,3]
            v = Materials[i_mat,4]
            y1 = Nodes[node_1,3]
            y2 = Nodes[node_2,3]
            h = abs(y2-y1)
            C = E*(1-v)/((1+v)*(1-2*v))*[1        v/(1-v)  0
                                         v/(1-v)  1        0
                                         0        0        (1-2*v)/(2*(1-v))]
            Me   =   kron(Hm,[rho*h 0; 0 rho*h])
            K1e  =   kron(Hm,[h*C[1,1] 0; 0 h*C[3,3]])
            K2e  = - kron(H2,[0 C[1,2]; C[3,3] 0]) + kron(H2',[0 C[3,3]; C[2,1] 0])
            K3e  =   kron(H3,[C[3,3]/h 0; 0 C[2,2]/h])
            DOFs =   collect(list_DOF[i_e,2]:list_DOF[i_e,3])
            K1g[DOFs,DOFs] = K1g[DOFs,DOFs] + K1e
            K2g[DOFs,DOFs] = K2g[DOFs,DOFs] + K2e
            K3g[DOFs,DOFs] = K3g[DOFs,DOFs] + K3e
            Mg[DOFs,DOFs]  = Mg[DOFs,DOFs]  + Me
        elseif Elements[i_e,2] == "2D_QuadTriangle"
            node_1 = Elements[i_e,5][1]
            node_2 = Elements[i_e,5][2]
            node_3 = Elements[i_e,5][3]
            node_4 = Elements[i_e,5][4]
            node_5 = Elements[i_e,5][5]
            node_6 = Elements[i_e,5][6]
            x = zeros(6)
            y = zeros(6)
            x[1] = Nodes[node_1,2]
            y[1] = Nodes[node_1,3]
            x[2] = Nodes[node_2,2]
            y[2] = Nodes[node_2,3]
            x[3] = Nodes[node_3,2]
            y[3] = Nodes[node_3,3]
            x[4] = Nodes[node_4,2]
            y[4] = Nodes[node_4,3]
            x[5] = Nodes[node_5,2]
            y[5] = Nodes[node_5,3]
            x[6] = Nodes[node_6,2]
            y[6] = Nodes[node_6,3]
            i_mat = Elements[i_e,3]
            E = Materials[i_mat,2]
            ρ = Materials[i_mat,3]
            ν = Materials[i_mat,4]

        end
        println("\nElement ", i_e, ":")
        println(node_1, " --> ", node_2, ": E = ", E, ", rho = ", rho, ", v = ", v, ", h = ", h)
        println("Evaluated without error.\n")
    end

    return K1g, K2g, K3g, Mg

end
