module TreeDynamics

include("math.jl")
include("tree.jl")
include("dynamics.jl")
include("plot.jl")

export TreeLinkage, TreeLimb, DynamicTree, make_tree!, collect_tree_nodes!, randtree!, all_leaf_paths, all_leaf_lines, ascend!, cache_end_position!
export cross, uncross, randr, r_euler1, r_euler2, r_euler3, isSO3, isso3, residualSO3, residualso3, axisangle
export dynamics!, proj_so3!, integrate_tree
export time_slide, plot_static_recursive_tree, do_recording

end
