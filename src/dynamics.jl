using LinearAlgebra
using ProgressMeter


function dynamics!(x::DynamicTree{<:Real}, t)
    M_c_root = zeros(3)
    M_k_root = zeros(3)
    M_c_leaf = zeros(3)
    M_k_leaf = zeros(3)
    
    relative_orientation = I(3)
    
    if !isnothing(AbstractTrees.parent(x))
        p = AbstractTrees.parent(x)
        # compute root-side moments
        relative_orientation = AbstractTrees.parent(x).limb.orientation'*x.limb.orientation
        
        relative_base_orientation = AbstractTrees.parent(x).limb.base_orientation'*x.limb.base_orientation
        
        # M_c_root = x.limb.damping*(relative_orientation'*AbstractTrees.parent(x).limb.rate - x.limb.rate)
        
        # this one seems good for damping:
        # M_c_root = x.limb.damping*(p.limb.rate - x.limb.rate)
        # M_c_root = x.limb.damping*(x.limb.orientation*p.limb.orientation'*p.limb.rate - x.limb.rate)
        
        # M_c_root = -x.limb.damping*x.limb.rate
        
        # M_k_root = -x.limb.stiffness*cross(relative_orientation*x.limb.direction, relative_base_orientation*x.limb.direction)
        # M_k_root = x.limb.stiffness*cross(x.limb.base_orientation*x.limb.direction, x.limb.orientation*x.limb.direction)
        # k_root_scale = p.limb.stiffness*norm(cross(p.limb.base_orientation'*p.limb.orientation*p.limb.direction, x.limb.base_orientation'*x.limb.orientation*x.limb.direction))
        # M_k_root = x.limb.stiffness*x.limb.orientation'*cross(cross(p.limb.orientation*p.limb.direction, p.limb.base_orientation*p.limb.direction), cross(x.limb.orientation*x.limb.direction, x.limb.base_orientation*x.limb.direction))
        M_k_root = x.limb.stiffness*x.limb.orientation'*(cross(p.limb.orientation*p.limb.direction, p.limb.base_orientation*p.limb.direction) .+ cross(x.limb.orientation*x.limb.direction, x.limb.base_orientation*x.limb.direction))
    end
    
    M_c_root = -x.limb.damping*x.limb.rate
    
    if !isempty(AbstractTrees.children(x))
        # compute leaf-side moments
        for tc in AbstractTrees.children(x)
            # first, cache all the child moments on this limb. Later will recurse.
            relative_orientation = tc.limb.orientation'*x.limb.orientation
            # relative_orientation = x.limb.orientation'*tc.limb.orientation
            
            relative_base_orientation = tc.limb.base_orientation'*x.limb.base_orientation
            # relative_base_orientation = x.limb.base_orientation'*tc.limb.base_orientation
            
            # M_c_leaf .+= x.limb.damping*(relative_orientation'*tc.limb.rate - x.limb.rate)
            
            # this one seems good for damping:
            # M_c_leaf .+= x.limb.damping*(tc.limb.rate - x.limb.rate)
            M_c_leaf .+= tc.limb.damping*(x.limb.orientation'*tc.limb.orientation*tc.limb.rate - x.limb.rate)
            # c_leaf_scale = norm(tc.limb.orientation*tc.limb.rate .- x.limb.orientation*x.limb.rate)
            
            # # if isnan(c_leaf_scale)
            # #     println(tc.limb.rate)
            # #     # println(tc.limb.orientation'*tc.limb.rate .- x.limb.orientation'*x.limb.rate)
            # # end
            # c_leaf_ax = zeros(3)
            # if norm(x.limb.rate) > 1e-9
            #     c_leaf_ax = x.limb.rate/norm(x.limb.rate)
            # end
            
            # M_c_leaf .+= c_leaf_scale * c_leaf_ax

            
            # M_c_leaf .+= x.limb.orientation*tc.limb.orientation'*tc.limb.damping*tc.limb.rate
            
            # M_k_leaf .+= -x.limb.stiffness*cross(relative_orientation'*x.limb.direction, relative_base_orientation'*x.limb.direction) 
            # M_k_leaf .+= x.limb.stiffness*cross(relative_orientation*x.limb.direction, relative_base_orientation*x.limb.direction) 
            # M_k_leaf .+= x.limb.stiffness*cross(x.limb.base_orientation*x.limb.direction, x.limb.orientation*x.limb.direction)
            # M_k_leaf = cross(tc.limb.base_orientation'*tc.limb.orientation*tc.limb.direction, x.limb.base_orientation'*x.limb.orientation*x.limb.direction)
            # M_k_leaf .+= tc.limb.stiffness*x.limb.orientation'*cross(cross(tc.limb.orientation*tc.limb.direction, tc.limb.base_orientation*tc.limb.direction), cross(x.limb.orientation*x.limb.direction, x.limb.base_orientation*x.limb.direction))
            M_k_leaf .+= tc.limb.stiffness*x.limb.orientation'*(cross(tc.limb.orientation*tc.limb.direction, tc.limb.base_orientation*tc.limb.direction) .+ cross(x.limb.orientation*x.limb.direction, x.limb.base_orientation*x.limb.direction))
        end
    end
    # println(norm(M_k_root), " ", norm(M_k_leaf))
    # println(norm(M_k_root), " ", norm(M_k_leaf), " ", norm(M_c_root), " ", norm(M_c_leaf))
    x.pocket["moment_root"] = M_c_root + M_k_root
    x.pocket["moment_leaf"] = M_c_leaf + M_k_leaf
    
    # at this point we can compute the forward step on this node. But want to do that in a separate pass to preserve the current attitude for descendent nodes to use.
    # Instead, just cache the dynamics here:
    
    x.pocket["deriv_orientation"] = -cross(x.limb.rate)*x.limb.orientation
    x.pocket["deriv_rate"] = pinv(x.limb.inertia)*(x.pocket["moment_root"] + x.pocket["moment_leaf"] - cross(x.limb.rate)*x.limb.inertia*x.limb.rate)
    
    for tc in AbstractTrees.children(x)
        dynamics!(tc,t)
    end
end

function step_tree_euler!(x::DynamicTree{<:Real}, dt, t)
    # note: 
    # I can advance a dcm C_BA by one timestep like this:
    #   C_BA[k+1] = exp(cross(omega)*dt)*C_BA[k]
    # 
    # The analytic timestepping of the angular acceleration equation (Euler) is a little more involved, I'm not sure I have the Jacobian correct.
    
    # new_orientation = exp(cross(x.limb.rate)*dt)*x.limb.orientation
    orientation_discrete_jacobian = exp(cross(x.limb.rate)*dt)
    new_orientation = x.limb.orientation*orientation_discrete_jacobian
    new_rate = x.limb.rate + x.pocket["deriv_rate"]*dt
    
    x.limb.orientation = new_orientation
    x.limb.rate = new_rate
    
    for tc in AbstractTrees.children(x)
        step_tree_euler!(tc, dt, t)
    end
end

function step_tree_rk4!(x::DynamicTree{<:Real}, dt, t)
    
end

function integrate_tree(x::DynamicTree{<:Real}, time::Vector{<:Real})
    origin = Point3f(0)
    nodes1 = [origin]
    collect_tree_nodes!(nodes1, x, [])
    
    record = [x]
    
    nodes = Matrix{Point3f}(undef, length(nodes1), length(time))
    orientations = zeros(3,3,length(nodes1) - 1, length(time))
    rates = zeros(3, length(nodes1) - 1, length(time))
    nodes[:,1] = nodes1
    @showprogress desc = "Integrating...\t" for (k,t) in enumerate(time)
        if k == 1
            continue
        end
        # to retain the whole thing, if desired:
        push!(record, deepcopy(x))
        # to retain only the nodes for plotting:
        #nodes = [Point3f(0.0)]
        #collect_tree_nodes(nodes x))
        
        # record the tree nodes
        nodesk = [origin]
        collect_tree_nodes!(nodesk, x, [])
        
        nodes[:,k] = nodesk
        
        # also record the DCMs and angular rates for each node in the tree
        for (i,node) in enumerate(AbstractTrees.PreOrderDFS(x))
            orientations[:,:,i,k] = node.limb.orientation
            rates[:,i,k] = node.limb.rate
        end
        
        dynamics!(x, t)     # to cache moments and derivatives jacobians
        
        # push the dynamics forward by a step:
        step_tree_euler!(x, time[k] - time[k - 1], t)
    end
    return record, nodes, orientations, rates
end

function dynamics!(dx::Vector{<:Real}, x::Vector{<:Real}, p::TreeLinkage{<:Real}, t)
    # function dynamics!(dx, x, p::TreeLinkage{<:Real}, t)
    # state: [w; C_Oi]
    # println("step ", t)
    z = [0; 0; 1]
    # this line should happen outside the dynamics function?

    w_i = reshape(x[1:3*p.n], 3, p.n)
    C_Oi = reshape(x[1+3*p.n:end], 3, 3, p.n)
    # is the reshaping happening in the wrong order?

    # re-orthonormalize dcm
    for k in 1:p.n
        (U, S, V) = svd(C_Oi[:, :, k])
        C_Oi[:, :, k] = U * diagm([1; 1; det(U) * det(V)]) * V'
    end

    # println("\tfirst dcm: ", C_Oi[:, :, 1])

    p.a_i[:, 1] = C_Oi[:, :, 1] * z
    p.a_i[:, 2] = C_Oi[:, :, 2] * z

    dC_Oi = zeros(3, 3, p.n)
    dw_i = zeros(3, p.n)
    for i in 2:p.n
        # moment from joint stiffness
        M_k = zeros(3, 1)
        M_c = zeros(3, 1)

        if i == p.n
            # last element
            M_k = p.k[i] * cross(p.a_i[:, i-1], p.a_i[:, i])
            # M_c = p.c[i] * (w_i[:, i-1] - w_i[:,i])
            M_c = -p.c[i] * w_i[:,i]
        else
            # other elements

            # this means i-1, i, i+1 have been updated by the time we get to this
            p.a_i[:, i+1] = C_Oi[:, :, i+1] * z

            M_k = p.k[i] * cross(p.a_i[:, i-1], p.a_i[:, i]) + p.k[i+1] * cross(p.a_i[:, i+1], p.a_i[:, i])
            # M_c = p.c[i] * (w_i[:, i-1] - w_i[:,i]) + p.c[i+1] * (w_i[:, i+1] - w_i[:,i])
            M_c = -p.c[i] * w_i[:,i] + p.c[i+1] * w_i[:, i+1]

        end

        # println(M_c)
        # angular acceleration
        dw_i[:, i] = pinv(p.I_i[:, :, i]) * (M_k + M_c - cross(w_i[:, i], p.I_i[:, :, i] * w_i[:, i]))

        # angle (attitude)
        dC_Oi[:, :, i] = -cross(w_i[:, i]) * C_Oi[:, :, i]

    end

    dx[:] = [vec(dw_i); vec(dC_Oi)]
end


# project the attitude states back to SO(3):
function proj_so3!(resid, x, p)
    C_Oi = reshape(x[1+3*p.n:end], 3, 3, p.n)
    resid[1:p.n] .= 0
    for i in 1:p.n
        bigI = 1 + 3 * p.n + 9 * (i - 1)
        # first DCM entry in resid gets residual, per DCM:
        # resid[bigI] = (det(C_Oi[:, :, i]) - 1)^2
        # all DCM entries in resid gets residual, per DCM:
        resid[bigI:bigI+8] .= (det(C_Oi[:, :, i]) - 1)^2 + (tr(C_Oi[:, :, i]' * C_Oi[:, :, i]) - 3)^2
    end
end
