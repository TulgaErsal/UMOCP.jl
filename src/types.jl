using JuMP
using DataFrames
using Parameters
using MadNLP
# These functions are required for NLOptMPC.jl and PrettyPlots.jl (resultsDir!)
# export  State,
#         Control,
#         Constraint,
#         Results,
#         Settings,
#         _Ipopt_defaults,
#         _Ipopt_MPC,
#         simulationModes,
#         MPC

# IPOPT Settings
# :print_level      : Print level
# :constr_viol_tol  : Absolute tolerance on the constraint violation.
# :max_iter         : Maximum number of iterations.
# :max_cpu_time     : A limit on CPU seconds that Ipopt can use to solve one problem
 _Ipopt_defaults=Dict{Symbol,Any}(
    :print_level                =>0,
    :warm_start_init_point      =>"yes",
    :tol                        =>1e-8,
    :max_iter                   =>3000,
    :max_cpu_time               =>1e6,
    :dual_inf_tol               =>1.,
    :constr_viol_tol            =>0.0001,
    :compl_inf_tol              =>0.0001,
    :acceptable_tol             =>1e-6,
    :acceptable_constr_viol_tol =>0.01,
    :acceptable_dual_inf_tol    =>1e-10,
    :acceptable_compl_inf_tol   =>0.01,
    :acceptable_obj_change_tol  =>1e20,
    :diverging_iterates_tol     =>1e20
)

_MadNLP_defaults=Dict{Symbol, Any}(
    # :linear_solver              => MadNLPMa27,
    :blas_num_threads           => 1,
    :print_level                => MadNLP.WARN,
    :tol                        => 1e-8,
    :acceptable_tol             => 1e-6,
    :max_iter                   => 3000,
    :max_wall_time              => 1e6,
    :jacobian_constant          => false,
    :hessian_constant           => false,
    :diverging_iterates_tol     => 1e20
)
#
# # IPOPT Settings for Model Predictive Control
 _Ipopt_MPC = Dict{Symbol,Any}(
    :print_level                =>0,
    :warm_start_init_point      =>"yes",
    :tol                        =>5e-1,
    :max_iter                   =>500,
    :max_cpu_time               =>0.47,
    :dual_inf_tol               =>5.,
    :constr_viol_tol            =>1e-1,
    :compl_inf_tol              =>1e-1,
    :acceptable_tol             =>1e-2,
    :acceptable_constr_viol_tol =>0.01,
    :acceptable_dual_inf_tol    =>1e10,
    :acceptable_compl_inf_tol   =>0.01,
    :acceptable_obj_change_tol  =>1e20,
    :diverging_iterates_tol     =>1e20
)
#
# # Simulation Modes
#
# ################################
# # Optimal Control Common Types #
# ################################
#
# # Control
@with_kw mutable struct Control
    name::Vector{Symbol}                = Vector{Symbol}[]          #
    description::Vector{AbstractString} = Vector{AbstractString}[]  #
    num::Int                            = Int(0)                    #
    pts::Int                            = Int(0)                    #
end
#
# # State
@with_kw mutable struct State
    name::Vector{Symbol}                 = Vector{Symbol}()          #
    description::Vector{AbstractString}  = Vector{AbstractString}()  #
    num::Int                             = 0                         #
    pts::Int                             = 0                         #
    model::Any                    = Any           #
end

# # Constraint
@with_kw mutable struct Constraint{ T <: Number }
    name::Vector{Symbol}                    = Vector{Symbol}()                  #
    handle::Vector{Any}      = Vector{Any}()      #
    value::Vector{Any}                        = Vector{Any}()                       #
    nums::Vector{Any}                       = Vector{Any}()                     # range of indices in g(x)
end
#
# # Solver
@with_kw mutable struct Solver
    name::Symbol                = :Ipopt            #
    settings::Dict{Symbol,Any}  = _Ipopt_defaults   #
end
#
# # Optimal Control Problem (OCP) Results
@with_kw mutable struct OCPResults{ T <: Number }
    tctr::Vector{Any}                           = Vector{T}()                           # Time vector for control
    tst::Vector{Any}                            = Vector{T}()                           # Time vector for state
    x           = Matrix{Any}[]      # JuMP states (scaled)
    u           = Matrix{Any}[]      # JuMP controls (scaled)
    xUnscaled::Matrix{VariableRef} = Matrix{VariableRef}(undef,0,0) # references to actual JuMP states
    uUnscaled::Matrix{VariableRef} = Matrix{VariableRef}(undef,0,0) # references to actual JuMP controls
    X                        = Matrix{T}[]                   # States
    U                        =  Matrix{T}[]                 # Controls
    X0                              = Vector{T}[]                           # Initial states for OCP
    CS                               = []                           # Costates
    tpolyPts                         = []                           # Time sample points for polynomials  (NOTE these interpolated solutions were developed for calculating error, between them and a known Optimal solution)
    XpolyPts                = []                   # State evaluated using Lagrange/Linear polynomial
    CSpolyPts                = []               # Costate evaluated using Lagrange/Linear polynomial
    UpolyPts                 = []                  # Control evaluated using Lagrane/Linear polynomial
    AlltpolyPts                      = []                           # Time sample points for polynomials
    AllXpolyPts              = []                  # State evaluated using Lagrange/Linear polynomial
    AllCSpolyPts             = []                   # Costate evaluated using Lagrange/Linear polynomial
    AllUpolyPts              = []                   # Control evaluated using Lagrane/Linear polynomial
    tpts::Vector{T}                             = Vector{T}()                           # Vector time sample points
    Xpts                    = Vector{T}[]                   # Vector state sample points
    Upts                     = Vector{T}[]                  # Vector control sample points
    CSpts                    = Vector{T}[]                  # Vector costate sample points
    x0Con          = nothing          # Handle for initial state constraints
    x0sCon          = nothing           # ? Unsure what this is yet (slack variable constraints?)
    xfCon           =  nothing          # Handle for final state constraints
    xfsCon          = nothing     # ? Unsure what this is yet (slack variable constraints?)
    dynCon = nothing   # Dynamics constraints
    constraint::Constraint{T}                   = Constraint{T}()                       # Constraint handles and data
    evalNum::Int                                = 1                                     # Number of times optimization has been run
    iterNum                                     = []                                    # Mics. data, perhaps an iteration number for a higher level algorithm
    status::MOI.TerminationStatusCode           = MOI.OTHER_ERROR                               # Optimization status
    tSolve::T                                   = convert(T,0)                          # Solve time for optimization
    objVal::T                                   = convert(T,0)                          # Objective function value
    dfs::Vector{DataFrame}                      = Vector{DataFrame}()                   # Results in DataFrame for plotting
    dfsOpt::DataFrame                           = DataFrame()                           # Optimization information in DataFrame for plotting
    dfsCon::DataFrame                           = DataFrame()                           # Constraint data

end
#
# # Optimal Control Problem (OCP) Settings
@with_kw mutable struct OCPSettings{ T <: Number }
    solver::Solver              = Solver()          # solver information
    finalTimeDV::Bool           = false             #
    integrationMethod::Symbol   = :ts               #
    integrationScheme::Symbol   = :bkwEuler         #
    save::Bool                  = true              # bool for saving data
    reset::Bool                 = false             # bool for reseting data
    evalConstraints::Bool       = false             # bool for evaluating duals of the constraints
    tfMin::T                    = convert(T,   0.0) # minimum final time
    tfMax::T                    = convert(T, 400.0) # maximum final time # ? Should choose something else as default?
    tfOptimal::Union{Bool, T}   = false             # known optimal final time
    numInterpPts::Int           = 250               # number of points to sample polynomial running through collocation points
    cacheOnly::Bool             = false             # bool for only caching the results when using optimize!()
    linearInterpolation::Bool   = false             # bool for using linear interpolation even if integrationMethod ==:ps
    interpolationOn::Bool       = false             # bool to indicate if user wants solution interpolated for them
    x0slackVariables::Bool      = false             #
    xFslackVariables::Bool      = false             #
    InternalLogging::Bool       = false             #
end
#
# # Optimal Control Problem (OCP) Flags
@with_kw mutable struct OCPFlags
  defined::Bool = false # a bool to tell if define() has been called
end
#
# # Optimal Control Problem (OCP)
@with_kw mutable struct OCP{T <: Number}

    # General properties
    state::State                            = State()                               # state data
    control::Control                        = Control()                             # control data
    tf             = Any                         # final time
    t0::JuMP.NonlinearParameter                      = @NLparameter(JuMP.Model(), x == 0)    # initial time # TODO: consider getting rid of this or replacing it with `n.mpc.v.t0Param`
    tV                           = Any                          # vector for use with time varying constraints

    # Boundary conditions
    X0::Vector{T}                           = Vector{T}()                           # initial state conditions
    X0_tol::Vector{T}                       = Vector{T}()                           # initial state tolerance
    x0s::Vector{JuMP.VariableRef}             = Vector{JuMP.VariableRef}()              # initial state variables
    XF::Vector{T}                           = Vector{T}()                           # final state conditions
    XF_tol::Vector{T}                       = Vector{T}()                           # final state tolerance
    xFs::Vector{JuMP.VariableRef}             = Vector{JuMP.VariableRef}()              # final state variables

    # Constant bounds on state variables
    XL::Vector{T}                           = Vector{T}()                           # Constant lower bound on state variables
    XU::Vector{T}                           = Vector{T}()                           # Constant upper bound on state variables

    # Variables for linear bounds on state variables
    mXL::Vector{Bool}                       = Vector{Bool}()                        # slope on XL -> time always starts at zero
    mXU::Vector{Bool}                       = Vector{Bool}()                        # slope on XU -> time always starts at zero
    XL_var = Any[]     # time varying lower bounds on states # ! NOTE: not used currently - was JuMP.Variable
    XU_var = Any[]     # time varying upper bounds on states # ! NOTE: not used currently - was JuMP.Variable

    # Constant bounds on control variables
    CL::Vector{T}                           = Vector{T}()                           # Constant lower bound on control variables
    CU::Vector{T}                           = Vector{T}()                           # Constant upper bound on control variables

    # Pseudospectral method data
    Nck::Vector{Int}                        = Vector{Int}()                         # number of collocation points per interval
    Nck_cum::Vector{Int}                    = Vector{Int}()                         # cumulative number of points per interval
    Nck_full::Vector{Int}                   = Vector{Int}()                         # [0;cumsum(n.ocp.Nck+1)]
    Ni::Int                                 = Int(0)                                # number of intervals
    tau::Array{Array{T,1},1}                = Array{Array{T,1},1}()                 # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
    ts::Array{Array{T,1},1} = Array{Array{T,1},1}()                                 # time scaled based off of tau
    w::Array{Array{T,1},1}  = Array{Array{T,1},1}()                                 # weights
    ws::Array{Array{T,1},1} = Array{Array{T,1},1}()                                 # scaled weights
    DMatrix::Vector{Matrix{T}}              = Vector{Matrix{T}}()                   # differention matrix
    IMatrix::Vector{Matrix{T}}              = Vector{Matrix{T}}()                   # integration matrix

    # tm method data
    N::Int                              = 0                                     # number of points in discretization
    dt::Array{Any,1}                    = Array{Any,1}()       # array of dts

    mdl::JuMP.Model                     = JuMP.Model()                          # JuMP model
    params                              = Any[]                                 # parameters for the models
    DXexpr                              = Any[]                                 # ? DX expression
    NLcon                               = Any[]                                 # ! NOTE: not used currently

    # Scaling factors
    XS::Vector{T}                       = Vector{T}()                           # scaling factors on states
    CS::Vector{T}                       = Vector{T}()                           # scaling factors on controls
end

#
@with_kw mutable struct Settings{ T <: Number }
    ocp::OCPSettings{T} = OCPSettings{T}()  #
end
#
# # Combined Results
@with_kw mutable struct Results{T <: Number}
    ocp::OCPResults{T}          = OCPResults{T}()           #
    resultsDir::AbstractString  = joinpath(pwd(),"results") # string that defines results folder
    mainDir::AbstractString     = pwd()                     # string that defines main folder
end
#
# # Combined Flags
@with_kw mutable struct Flags
    ocp::OCPFlags   = OCPFlags()    #
end
#
# # Nonlinear Optimal Control (NLOpt) Problem
@with_kw mutable struct NLOpt{T <: Number}
    # major data structs
    ocp::OCP{T}     = OCP{T}()      #
    s::Settings     = Settings{T}() #
    r::Results{T}   = Results{T}()  #
    f::Flags        = Flags()       #
end
#
# # Default Number subtype is Float64
NLOpt_realization() = NLOpt{Float64}()
