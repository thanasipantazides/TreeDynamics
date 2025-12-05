using TreeDynamics
using GLMakie
using LinearAlgebra
import AbstractTrees

function main()
    GLMakie.activate!()
    
    fig = GLMakie.Figure(size=(400,400))
    fig = Figure(size=(600,600))
    # display(fig)
    layout = GridLayout(fig[1, 1], tellheight = false)
    # ax2d = Axis(fig[2,1])
    
    display(fig, title="TreeDynamics")
    
    inertia = [
        1 0 0;
        0 1 0;
        0 0 4] * 1e-2
    t_root = DynamicTree(TreeLimb(
        1.0,        # length
        [0.0, 0.0, 1.0],    # direction
        2.0e-2,        # stiffness
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
    randtree!(7,3,t_root, length=1)
    # randtree!(5,5,t_root)
    
    for leaf in AbstractTrees.PreOrderDFS(t_root)
        if isempty(AbstractTrees.children(leaf))
            leaf.pocket["perturb"] = 1
        end
    end
    println("Animating...")
    dt = 0.05
    idle_animate(layout, t_root, dt)
    
    
    # phase = 0
    # x = 0:0.1:2*pi
    # pt = Observable(Point2f[(0.0,0.0) for k in 1:length(x)])
    
    # lines!(ax2d, pt)
    # limits!(ax2d, 0,2*pi, -1, 1)
    
    # fps = 60
    
    # while true
    #     pt[] = [Point2f(xi, cos(xi + phase)) for xi in x]
    #     phase = (phase + 0.1)
    #     notify(pt)
    #     sleep(1/fps)
    # end
end