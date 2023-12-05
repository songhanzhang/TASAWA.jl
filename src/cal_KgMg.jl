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
            node_1 = Int(Elements[i_e,5][1])
            node_2 = Int(Elements[i_e,5][2])
            node_3 = Int(Elements[i_e,5][3])
            node_4 = Int(Elements[i_e,5][4])
            node_5 = Int(Elements[i_e,5][5])
            node_6 = Int(Elements[i_e,5][6])
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
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            ν = Materials[i_mat,2][3]
            (K1e, K2e, K3e, Me) = cal_KeMe_QuadraticTriangular(x,y,E,ρ,ν)
            DOF_1  = Int(findall(isequal(node_1+0.1),list_DOF[:,2])[1])
            DOF_2  = Int(findall(isequal(node_1+0.2),list_DOF[:,2])[1])
            DOF_3  = Int(findall(isequal(node_1+0.3),list_DOF[:,2])[1])
            DOF_4  = Int(findall(isequal(node_2+0.1),list_DOF[:,2])[1])
            DOF_5  = Int(findall(isequal(node_2+0.2),list_DOF[:,2])[1])
            DOF_6  = Int(findall(isequal(node_2+0.3),list_DOF[:,2])[1])
            DOF_7  = Int(findall(isequal(node_3+0.1),list_DOF[:,2])[1])
            DOF_8  = Int(findall(isequal(node_3+0.2),list_DOF[:,2])[1])
            DOF_9  = Int(findall(isequal(node_3+0.3),list_DOF[:,2])[1])
            DOF_10 = Int(findall(isequal(node_4+0.1),list_DOF[:,2])[1])
            DOF_11 = Int(findall(isequal(node_4+0.2),list_DOF[:,2])[1])
            DOF_12 = Int(findall(isequal(node_4+0.3),list_DOF[:,2])[1])
            DOF_13 = Int(findall(isequal(node_5+0.1),list_DOF[:,2])[1])
            DOF_14 = Int(findall(isequal(node_5+0.2),list_DOF[:,2])[1])
            DOF_15 = Int(findall(isequal(node_5+0.3),list_DOF[:,2])[1])
            DOF_16 = Int(findall(isequal(node_6+0.1),list_DOF[:,2])[1])
            DOF_17 = Int(findall(isequal(node_6+0.2),list_DOF[:,2])[1])
            DOF_18 = Int(findall(isequal(node_6+0.3),list_DOF[:,2])[1])
            DOFs = [DOF_1,DOF_2,DOF_3,DOF_4,DOF_5,DOF_6,DOF_7,DOF_8,DOF_9,DOF_10,DOF_11,DOF_12,DOF_13,DOF_14,DOF_15,DOF_16,DOF_17,DOF_18]
            K1g[DOFs, DOFs] += K1e
            K2g[DOFs, DOFs] += K2e
            K3g[DOFs, DOFs] += K3e
            Mg[DOFs, DOFs] += Me
            println("\nElement ", i_e, ":")
            println("[", node_1, ",", node_2, ",", node_3, ",", node_4, ",", node_5, ",", node_6, "]: E = ", E, ", rho = ", ρ, ", v = ", ν)
        end
        println("\nEvaluated without error.\n")
    end

    return K1g, K2g, K3g, Mg

end
