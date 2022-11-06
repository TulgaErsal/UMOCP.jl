isdefined(Base, :__precompile__) && __precompile__()

module UMOCP

using JuMP
import JuMP: value,
             @NLexpression,
             @NLobjective,
             @NLparameter,
             @NLconstraint,
             Model
export @NLexpression,
       @NLobjective,
       @NLparameter,
       @NLconstraint,
       value,
       Model


using FastGaussQuadrature
using DataFrames
using Parameters
using MadNLP
using Ipopt

# Write your package code here.
include("math.jl")
export lagrange_basis_poly
export interpolate_lagrange
export scale_w
export scale_tau
export lagrange_basis_poly!
export lagrange_basis_poly
export polyDiff
export linearSpline

include("types.jl")

include("utils.jl")
export intervals
export initState
export initControl
export states!
export controls!
export defineTolerances!
export create_tV!
export initConstraint!
export newConstraint!
export integrate!
export NLoptimize!
export evalConstraints!
export postProcess!
export opt2dfs!
export dvs2dfs
export try_import

include("setup.jl")
export define
export defineSolver!
export OCPdef!
export configure!

include("ps.jl")
export DMatrix!
export createIntervals!

include("diffeq.jl")
export dynamics!
export constraints!
export DiffEq
export addCon
export NLExpr
export NLCon



include("parameter.jl")
include("ThreeDOF_Bicycle.jl")
export ThreeDOFBicycle_expr
export ThreeDOF_clay_nn_expr

end
