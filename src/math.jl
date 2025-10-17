using LinearAlgebra
using StaticArrays


so3tol = 1e-4
    
function isSO3(X::AbstractMatrix)::Bool
    if (det(X) - 1 < so3tol) && (abs(sum(X'*X .- diagm([1,1,1]))) < so3tol)
        return true
    end
    return false
end

function residualSO3(X::AbstractMatrix)
    return (det(X) - 1) + (sum(X'*X .- diagm([1,1,1])))
end

function residualso3(X::AbstractMatrix)
    return sum(X' + X) 
end

function isso3(X::AbstractMatrix)::Bool
    if abs(sum(X' + X)) < so3tol
        return true
    end
    return false
end

function r_euler3(ang::Real)::SMatrix{3,3}{<:Real}
    return [cos(ang) -sin(ang) 0;
        sin(ang) cos(ang) 0;
        0 0 1]
end
function r_euler2(ang::Real)::SMatrix{3,3}{<:Real}
    return [cos(ang) 0 sin(ang);
        -sin(ang) 0 cos(ang);
        0 1 0]
end
function r_euler1(ang::Real)::SMatrix{3,3}{<:Real}
    return [1 0 0;
        0 cos(ang) -sin(ang);
        0 sin(ang) cos(ang)]
end
import LinearAlgebra.cross

function LinearAlgebra.cross(x::AbstractVector)
    if length(x) > 3
        throw(Exception("too long!"))
    end

    return [0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0]
end

function uncross(X::AbstractMatrix)
    return [-X[2,3]; X[1,3]; -X[1,2]]
end

function randr()::Matrix{<:Real}
    ax = rand(3) .- 1
    ax = ax / norm(ax)
    cang = 2*rand() - 1

    return I * cang + (1 - cang) * ax * ax' + cross(ax) * sqrt(1 - cang^2)
end