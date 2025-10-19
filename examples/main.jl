using TreeDynamics
using GLMakie
using LinearAlgebra
import AbstractTrees
import Base: show

function Base.show(io::IO, n::DynamicTree)
   str_parent = isnothing(n.parent) ? "no parent" : "with parent"
   str_children = "$(length(AbstractTrees.children(n))) children"
   print(io, "DynamicTree($str_parent, $(n.limb.stiffness), $(n.limb.damping), $str_children)")
end

function main()
    inertia = [
        10 0 0;
        0 10 0;
        0 0 40] * 1e-3
    t_root = DynamicTree(TreeLimb(
        1.0,        # length
        [0.0, 0.0, 1.0],    # direction
        1.0e-2,        # stiffness
        1e-2,       # damping
        inertia,    # inertia
        # randr(),    # base orientation
        # r_euler1(pi/6),
        I(3),
        I(3),       # orientation (dynamic)
        [0.0, 0.0, 0.01] # angular rate
    ), [], nothing, nothing)
    
    println("Building tree...")
    # make_tree!(5,3,t_root)
    randtree!(4,3,t_root)
    t_root_init = deepcopy(t_root)
    
    startnodes = [Point3f(0.0)]
    color_root = [RGBAf(0.0,0.0,0.0,1.0)]
    collect_tree_nodes!(startnodes, t_root, color_root)
    println("\tlength of startnodes: ", length(startnodes))
    treesize = AbstractTrees.treesize(t_root)
    println("\tnumber of nodes: ", treesize)
    
    println("Integrating...")
    step = 0.01
    time = Vector(0:step:10)
    trees, nodes, orientations, rates = integrate_tree(t_root, time)
        
    println("Checking health...")
    
    bSO3 = zeros(treesize, length(time))
    ratesnorms = zeros(treesize, length(time))
    for k in 2:length(time)
        for i in 1:treesize
            bSO3[i,k] = residualSO3(orientations[:,:,i,k])
            ratesnorms[i,k] = norm(rates[:, i, k])
        end
    end
        
    
    println("Plotting...")
    
    GLMakie.activate!()
    fig = Figure(size=(600,600))
    # display(fig)
    layout = GridLayout(fig[1, 1], tellheight = false)
    # al = AmbientLight(RGBf(0.4,0.4,0.4))
    layout3d = GridLayout(layout[1,1], tellheight=false)
    ax3d, slider = time_slide(layout3d, trees, time)
    # do_recording(ax3d, slider, time, "doc/assets/tree_movie.mp4")
    display(fig)
    
    layout2d = GridLayout(layout[1,2])
    
    # ls,cs = all_leaf_lines(t_root_init)
    # nl = length(ls)
    # for i in eachindex(ls)
    #     lines!(ax3d, ls[i], color=cs[i])
    # end  
    # # lines!(ax3d, startnodes, color=RGBAf(0.4,0.73,0.0,1.0), linewidth=0.5)
    # # lines!(ax3d, startnodes, color=color_root, linewidth=0.5)
    # ##########################################
    # return t_root
    
    axSO3 =  Axis(layout2d[1,1], xlabel="Time", ylabel=L"\mathrm{SO(3)-ness}")
    axrate = Axis(layout2d[2,1], xlabel="Time", ylabel=L"\mathbf{\omega}_B")
    
    linkxaxes!(axSO3, axrate)
    
    for i in 1:treesize
        lines!(axSO3, time, bSO3[i,:])
        # lines!(axrate, time, ratesnorms[i,:])
        lines!(axrate, time, rates[1,i,:], color=:black, linewidth=0.25)
        lines!(axrate, time, rates[2,i,:], color=:black, linewidth=0.25)
        lines!(axrate, time, rates[3,i,:], color=:black, linewidth=0.25)
    end
    
    return t_root
end