import AbstractTrees
using GLMakie

function idle_animate(layout::GLMakie.GridLayout, tree::DynamicTree, dt::Real)
    al = AmbientLight(RGBf(0.4,0.4,0.4))
    ax3d = LScene(
        layout[1,1],
        show_axis = false,
        scenekw = (
            # lights=[dl, al],
            lights = [al],
            # backgroundcolor=:black,
            clear=true
        ),
        tellheight = false
    )
    n_nodes = length(collect(AbstractTrees.PreOrderDFS(tree)))
    nobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:n_nodes])
    cobs = Observable([RGBAf(0.0,0.0,0.0,0.0) for i in 1:n_nodes])
    wobs = Observable([1.0 for i in 1:n_nodes])
    lines!(ax3d, nobs, color=cobs, linewidth=wobs)
    
    fps = 60
    originp = Point3f(0.0)
    originc = RGBAf(0,0,0,1)
    originw = 8.0
    t = 0
    while true # integrate the model forever
        t += dt
        dynamics!(tree, t)
        step_tree_euler!(tree, dt, t)
        
        nodesk = [originp]
        colsk = [originc]
        widthsk = [originw]
        collect_tree_nodes!(nodesk, tree, colsk, widthsk)
        nobs[] = nodesk
        cobs[] = colsk
        wobs[] = widthsk
        notify(nobs)
        notify(cobs)
        notify(wobs)
        
        # cache_end_position!(originp, tree)
        
        sleep(1/fps)
    end
end

function time_slide(layout::GLMakie.GridLayout, trees::Vector{DynamicTree{Float64}}, time::Vector{<:Real})
    al = AmbientLight(RGBf(0.4,0.4,0.4))
    ax3d = LScene(
        layout[1,1],
        show_axis = false,
        scenekw = (
            # lights=[dl, al],
            lights = [al],
            # backgroundcolor=:black,
            clear=true
        ),
        tellheight = false
    )
    slide = Slider(layout[2,1], range=eachindex(time), startvalue=1)
    
    on(events(ax3d.scene).scroll, priority=100) do (dx, dy)
        if ispressed(ax3d, Keyboard.left_super)
            next = Int(slide.value[] + 5*sign(dy))
            next = min(max(next, 1), length(time))
            # slide.value[] = next
            set_close_to!(slide, next)
            # notify(slide.value)
            Consume(true)
        end
    end
    
    n_nodes = length(collect(AbstractTrees.PreOrderDFS(trees[1])))
    
    nobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:n_nodes])
    cobs = Observable([RGBAf(0.0,0.0,0.0,0.0) for i in 1:n_nodes])
    
    xstaobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    xdynobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    ystaobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    ydynobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    zstaobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    zdynobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    
    momrootobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    momleafobs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:3*n_nodes])
    
    originp = Point3f(0.0)
    originc = RGBAf(0,0,0,1)
    on(slide.value) do k
        nodesk = [originp]
        colsk = [originc]
        collect_tree_nodes!(nodesk, trees[k], colsk)
        nobs[] = nodesk
        cobs[] = colsk
        notify(nobs)
        notify(cobs)
        
        cache_end_position!(originp, trees[k])
        
        i3 = 1
        xdynobs[] = [Point3f(0.0)]
        ydynobs[] = [Point3f(0.0)]
        zdynobs[] = [Point3f(0.0)]
        momrootobs[] = [Point3f(0.0)]
        momleafobs[] = [Point3f(0.0)]
        nodes = collect(AbstractTrees.PreOrderDFS(trees[k]))
        for i in 1:n_nodes
            # origin position
            dcm_origin = Point3f(nodes[i].pocket["end_position"] - nodes[i].limb.length/2*nodes[i].limb.orientation*nodes[i].limb.direction)
            # x-direction end position
            new_x = Point3f(dcm_origin + nodes[i].limb.length/4*nodes[i].limb.orientation[:,1])
            new_y = Point3f(dcm_origin + nodes[i].limb.length/4*nodes[i].limb.orientation[:,2])
            new_z = Point3f(dcm_origin + nodes[i].limb.length/4*nodes[i].limb.orientation[:,3])
            append!(xdynobs[], [dcm_origin, new_x, Point3f(NaN)])
            append!(ydynobs[], [dcm_origin, new_y, Point3f(NaN)])
            append!(zdynobs[], [dcm_origin, new_z, Point3f(NaN)])
            
            momleaf = nodes[i].pocket["moment_leaf"]
            momroot = nodes[i].pocket["moment_root"]
            new_momleaf = norm(momleaf) > 1e-9 ? Point3f(dcm_origin + 2*nodes[i].limb.length/4*nodes[i].limb.orientation*momleaf./norm(momleaf)) : dcm_origin
            new_momroot = norm(momroot) > 1e-9 ? Point3f(dcm_origin + 2*nodes[i].limb.length/4*nodes[i].limb.orientation*momroot./norm(momroot)) : dcm_origin
            append!(momrootobs[], [dcm_origin, new_momroot, Point3f(NaN)])
            append!(momleafobs[], [dcm_origin, new_momleaf, Point3f(NaN)])
            # i3 += 3
        end
        notify(xdynobs)
        notify(ydynobs)
        notify(zdynobs)
        notify(momrootobs)
        notify(momleafobs)
        # orientation:
        # How to do this without recursion? Would be extra handy if 
        # xstaobs = []
    end
    lines!(ax3d, nobs, color=cobs, linewidth=2)
    
    lines!(ax3d, xdynobs, color=:red, linewidth=0.5)
    lines!(ax3d, ydynobs, color=:green, linewidth=0.5)
    lines!(ax3d, zdynobs, color=:blue, linewidth=0.5)
    lines!(ax3d, momrootobs, color=:purple, linewidth=0.5)
    lines!(ax3d, momleafobs, color=:cyan, linewidth=0.5)
    
    return ax3d, slide
end

function do_recording(ax::GLMakie.LScene, slide::GLMakie.Slider, time::Vector{<:Real}, file::String)
    println("Saving video...")
    framerate = 30
    timestamps = range(1, length(time), step=2)
    record(ax.parent, file, timestamps; framerate = framerate) do tk
        set_close_to!(slide, tk)
    end
end

function time_slide(fig::GLMakie.Figure, sol, p::TreeLinkage)
    a = [0; 0; 1]

    C_Oih = reshape(sol[3*p.n+1:end, :], 3, 3, p.n, length(sol.t))

    # joint positions and member lengths
    r = zeros(3, p.n, length(sol.t))
    rn = zeros(p.n, length(sol.t))
    Cdet = zeros(p.n, length(sol.t))
    Ctr = zeros(p.n, length(sol.t))
    for t in 1:length(sol.t)
        for i = 1:p.n
            if i > 1
                r[:, i, t] = r[:, i-1, t] + p.l[i] * C_Oih[:, :, i, t] * a
                rn[i, t] = norm(r[:, i, t])
            end
            Cdet[i,t] = det(C_Oih[:,:,i,t])
            Ctr[i,t] = tr(C_Oih[:,:,i,t]'*C_Oih[:,:,i,t]) - 2
        end
    end

    layout = GridLayout(fig[1, 1], tellheight = false)
    al = AmbientLight(RGBf(0.4,0.4,0.4))
    layout3d = GridLayout(layout[1,1], tellheight=false)
    layout2d = GridLayout(layout[1,2])
    ax3d = LScene(
        layout3d[1,1],
        show_axis = true,
        scenekw = (
            # lights=[dl, al],
            lights = [al],
            # backgroundcolor=:black,
            clear=true
        ),
        tellheight = false
    )
    axCdet = Axis(layout2d[1,1],xlabel="time", ylabel=L"\det(C)")
    axCtr = Axis(layout2d[2,1],xlabel="time", ylabel=L"\mathrm{tr}(C^TC) - 2 ")

    println(eachindex(sol.t))
    slide = Slider(layout3d[2,1], range=eachindex(sol.t), startvalue=1)

    hlen = Int(floor(length(sol.t)/10))
    # hlen = 100
    robs = Observable(Point3f[(0.0,0.0,0.0) for i in 1:p.n])
    # rhist = Observable(Point3f[NaN for i in 1:p.n])
    rhist = Observable([Point3f(NaN)])
    chist = Observable([RGBAf(0.0,0.0,0.0,0.0)])

    on(slide.value) do k
        for i in 1:p.n
            robs[][i] = r[:,i,k]
        end

        rhist[] = [Point3f(NaN)]
        chist[] = [RGBAf(0.0,0.0,0.0,0.0)]
        startk = max(k - hlen, 1)
        L = length(startk:k)
        for j in startk:k
            append!(rhist[], Point3f[r[:,i,j] for i in 1:p.n])
            append!(chist[], [RGBAf(0.4*(p.n - i)/p.n,0.73*(j - startk)/L,0.0,Float32(j - startk)/L) for i in 1:p.n])
            push!(chist[], RGBAf(0.4,0.73,0.0,0.0))
            push!(rhist[], Point3f(NaN))
        end
        notify(chist)
        notify(rhist)
        notify(robs)
    end
    tslide = lift(slide.value) do v
        return sol.t[v]
    end

    scatterlines!(ax3d, robs, color=RGBAf(0.4,0.73,0.0,1.0))
    lines!(ax3d, rhist, color=chist, linewidth=0.25)

    lines!(ax3d, [Point3f(0.0), Point3f(0.0,0.0,sum(p.l))], color=:black, alpha=0.25, linewidth=0.1)

    [lines!(axCdet, sol.t, Cdet[i,:]) for i in 1:p.n]
    [lines!(axCtr, sol.t, Ctr[i,:]) for i in 1:p.n]
    vlines!(axCdet, tslide, color=RGBAf(0.4,0.73,0.0,1.0))
    vlines!(axCtr, tslide, color=RGBAf(0.4,0.73,0.0,1.0))
end

function plot_static_recursive_tree(ax::GLMakie.LScene, p::DynamicTree{<:Real})
    if isempty(AbstractTrees.children(p))
        return
    end
    p0 = Point3f(0.0) # an inital point to root the plot
    if isnothing(AbstractTrees.parent(p))
        # if the parent is nothing, this node is root. Cache its position for descendents to use.
        p.pocket["end_position"] = p0 + Point3f(p.limb.length*p.limb.base_orientation*p.limb.direction)
    else
        # if there is a parent, use their cached position later.
        p0 = p.pocket["end_position"]
    end
    
    for tc in AbstractTrees.children(p)
        # start of this child branch is p0 (the parent's root) + the parent's length, along their direction
        root = p0
        # end of this child branch is the start of it + its length along its direction
        tip = p0 .+ Point3f(tc.limb.length*tc.limb.base_orientation*tc.limb.direction)
        # draw it
        lines!(ax, [p0, tip], color=RGBAf(0.4,0.73,0.0,1.0))
        # cache this position for descendents
        tc.pocket["end_position"] = tip
        plot_static_recursive_tree(ax, tc)
    end
end