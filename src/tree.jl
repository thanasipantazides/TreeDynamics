using StaticArrays
using GLMakie
import AbstractTrees: children, nodevalue, parent, childtype, ParentLinks, StoredParents
# using AbstractTrees: children, nodevalue, parent, childtype, ParentLinks, StoredParents

mutable struct TreeLinkage{T<:Real}
    n::Int64            # static
    l::Vector{T}        # .
    k::Vector{T}        # .
    c::Vector{T}        # .
    I_i::Array{T,3}     # .
    a_i::Matrix{T}      # .
    # C_Oi::Array{T,3}    # dynamic
    # w_i::Matrix{T}      # dynamic
end

# a container for static and dynamic properties of individual limbs.
mutable struct TreeLimb{T<:Real}
    length::T                       # branch length
    direction::SVector{3}{T}        # branch direction in the local coordinate system
    stiffness::T                    # spring constant for this-to-parent joint
    damping::T                      # damping coefficient for this-to-parent joint
    inertia::SMatrix{3,3}{T}        # inertia matrix
    base_orientation::SMatrix{3,3}{T}   # static rotation transform w.r.t inertial.
    orientation::SMatrix{3,3}{T}    # rotation transform w.r.t. inertial
    rate::SVector{3}{T}             # angular rate, dynamic
end

function TreeLimb(
    length,
    direction::AbstractVector{<:Real},
    stiffness,
    damping,
    inertia::AbstractMatrix{<:Real},
    base_orientation::AbstractMatrix{<:Real},
    orientation::AbstractMatrix{<:Real},
    rate::AbstractVector{<:Real}
)
    TS = typeof(promote(length, stiffness, damping)[1])
    TA = Base.promote_eltype(direction, inertia, base_orientation, orientation, rate)
    TR = promote_type(TS, TA)
    return TreeLimb(
        convert(TR, length),
        convert(SVector{3}{TR}, direction),
        convert(TR, stiffness),
        convert(TR, damping),
        convert(SMatrix{3,3}{TR}, inertia),
        convert(SMatrix{3,3}{TR}, base_orientation),
        convert(SMatrix{3,3}{TR}, orientation),
        convert(SVector{3}{TR}, rate)
    )
end

mutable struct DynamicTree{T<:Real}
    limb::TreeLimb{T}
    children::Vector{DynamicTree{T}}
    parent::Union{Nothing,DynamicTree{T}}
    pocket::Dict{String, Any}
    function DynamicTree(limb::TreeLimb{T}, children::Any, parent=nothing, pocket::Union{Nothing,Dict{String, Any}}=nothing) where T<:Real
        if isempty(children)
            children = DynamicTree[]
        end
        if isnothing(pocket)
            pocket = Dict{String, Any}()
        end
        new{T}(limb, children, parent, pocket)
    end
end

# methods to implement the AbstractTree interface:
children(t::DynamicTree) = t.children
nodevalue(t::DynamicTree) = t.limb
childtype(t::DynamicTree) = TreeLimb
ParentLinks(::Type{<:DynamicTree{<:Real}}) = StoredParents()
parent(n::DynamicTree) = n.parent

# recursive construction of trees:
function make_tree!(n_depth::Int, n_breadth::Int, tree::Union{Nothing,DynamicTree{<:Real}}=nothing)
    if n_depth == 0
        return
    end
    
    for i in 1:n_breadth
        # do stuff to modify tc
        direction = tree.limb.direction
        stiffness = 2/3*tree.limb.stiffness
        damping = 2/3*tree.limb.damping
        
        length = 1/2*tree.limb.length
        
        rate = zeros(3)
        
        circang = 2*pi*(i - 1)/n_breadth
        skewang = pi/4
        transform = r_euler3(circang)*r_euler1(skewang)
        base_orientation = tree.limb.base_orientation*transform
        orientation = base_orientation
        inertia = 0.75*tree.limb.inertia
        
        tc = TreeLimb(
            length,
            direction,
            stiffness,
            damping,
            inertia,
            base_orientation,
            orientation,
            rate
        )
        
        nt = DynamicTree(tc, [], tree, nothing)
        
        push!(tree.children, nt)
        make_tree!(n_depth - 1, n_breadth, nt)
    end
end

function randtree!(n_depth::Int, n_breadth::Int, tree::Union{Nothing,DynamicTree{<:Real}}=nothing)
    if n_depth == 0
        return
    end
    
    for i in 1:n_breadth
        # do stuff to modify tc
        direction = tree.limb.direction
        stiffness = 1/3*tree.limb.stiffness
        damping = 1/3*tree.limb.damping
        length = 1/2*tree.limb.length
        
        rate = rand(3)*1e0
        # rate = zeros(3)
        
        circang = 2*pi*(i - 1)/n_breadth #+ rand()*pi/10
        skewang = pi/4 #+ rand()*pi/12
        transform = r_euler3(circang)*r_euler1(skewang)
        base_orientation = tree.limb.base_orientation*transform
        # mv_orientation = r_euler1(pi/6)
        # orientation = base_orientation*mv_orientation
        orientation = base_orientation
        # orientation = randr() 
        inertia = tree.limb.inertia*0.66
        
        tc = TreeLimb(
            length,
            direction,
            stiffness,
            damping,
            inertia,
            base_orientation,
            orientation,
            rate
        )
        
        nt = DynamicTree(tc, [], tree, nothing)
        
        push!(tree.children, nt)
        make_tree!(n_depth - 1, n_breadth, nt)
    end
end

function ascend!(tree::DynamicTree{<:Real}, path)
    push!(path, tree)
    if isnothing(AbstractTrees.parent(tree))
        return
    else
        ascend!(AbstractTrees.parent(tree), path)
    end
end

function all_leaf_paths(tree::DynamicTree{<:Real})
    nleaves = length(collect(AbstractTrees.Leaves(tree)))
    
    paths = []
    
    for l in collect(AbstractTrees.Leaves(tree))
        p = []
        ascend!(l, p)
        push!(paths, reverse!(p))
    end
    return paths
end

function all_leaf_lines(tree::DynamicTree{<:Real})
    paths = all_leaf_paths(tree)
    
    root = Point3f(0.0)
    
    lines = []
    colors = []
    for p in paths
        line = [root]
        color = [RGBAf(0.0,0.0,0.0,1.0)]
        for n in p
            point = Point3f(0.0)
            if isnothing(AbstractTrees.parent(n))
                point = root + Point3f(n.limb.length*n.limb.orientation*n.limb.direction)
                n.pocket["end_position_alt"] = point
                
                c = RGBAf(1.0,0.0,0.0,1.0)
            else
                elder = AbstractTrees.parent(n)
                
                point = n.parent.pocket["end_position_alt"] + Point3f(n.limb.length*n.limb.orientation*n.limb.direction)
                n.pocket["end_position_alt"] = point
                
                c = RGBAf(0.0,1.0,0.0,1.0)
            end
            push!(line, point)
            push!(color, c)
        end
        push!(lines, line)
        push!(colors, color)
    end
    return lines, colors
end

function cache_end_position!(root::Point3f, tree::DynamicTree{<:Real})
    for t in AbstractTrees.PreOrderDFS(tree)
        if isnothing(AbstractTrees.parent(t))
            t.pocket["end_position"] = root + Point3f(t.limb.length*t.limb.orientation*t.limb.direction)
        else
            t.pocket["end_position"] = t.parent.pocket["end_position"] + Point3f(t.limb.length*t.limb.orientation*t.limb.direction)
        end
    end
end

function collect_tree_nodes!(root::Vector{Point3f}, node::DynamicTree{<:Real}, color_root::Vector{<:Any})
    # trawl the tree and add points for quick plotting
    
    for t in AbstractTrees.PreOrderDFS(node)
        if isnothing(AbstractTrees.parent(t))
            # if this is the root node, draw again from the initial root element (root[1]):
            # t.pocket["end_position"] = root[1] + Point3f(t.limb.length*t.limb.base_orientation*t.limb.direction)
            t.pocket["end_position"] = root[1] + Point3f(t.limb.length*t.limb.orientation*t.limb.direction)
            
            push!(root, t.pocket["end_position"])
            push!(color_root, RGBAf(174/255, 91/255, 46/255, 1.0))
        else
            # if this has a parent, draw starting from the parent's end position:
            # t.pocket["end_position"] =  t.parent.pocket["end_position"] + Point3f(t.limb.length*t.limb.base_orientation*t.limb.direction)
            t.pocket["end_position"] = t.parent.pocket["end_position"] + Point3f(t.limb.length*t.limb.orientation*t.limb.direction)
            
            # push!(root, t.parent.pocket["end_position"])
            push!(root, t.parent.pocket["end_position"])
            push!(color_root, RGBAf(0.4,0.73,0.0,1.0))
            push!(root, t.pocket["end_position"])
            push!(color_root, RGBAf(0.4,0.73,0.0,1.0))
        end
        
        if isempty(AbstractTrees.children(t))
            push!(root, Point3f(NaN))
            push!(color_root, RGBAf(0.4,0.73,0.0,0.0))
        end
    end
end