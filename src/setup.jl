using JuMP
# using FastGaussQuadrature
# using Ipopt

function define(;
    numStates::Int = 0,
    numControls::Int = 0,
    X0 = fill(NaN,numStates),
    XF = fill(NaN,numStates),
    XL = fill(NaN,numStates),
    XU = fill(NaN,numStates),
    CL = fill(NaN,numControls),
    CU = fill(NaN,numControls)
)::NLOpt

    n = NLOptRealization()

    # validate input
    if numControls <= 0
        error("Controls ($numControls) must be > 0")
    end
    if numStates <= 0
        error("States ($numStates) must be > 0")
    end
    if length(X0) != numStates
        error("Length of X0 ($(length(X0))) must match number of states ($numStates)")
    end
    if length(XF) != numStates
        error("Length of XF ($(length(XF))) must match number of states ($numStates)")
    end
    if length(XL) != numStates
        error("Length of XL ($(length(XL))) must match number of states ($numStates)")
    end
    if length(XU) != numStates
        error("Length of XU ($(length(XU))) must match number of states ($numStates)")
    end
    if length(CL) != numControls
        error("Length of CL ($(length(CL))) must match number of controls ($numControls)")
    end
    if length(CU) != numControls
        error("Length of CU ($(length(CU))) must match number of controls ($numControls)")
    end

    n.ocp.state     = initState(numStates)
    n.ocp.control   = initControl(numControls)
    n.ocp.X0        = X0
    n.ocp.X0_tol    = fill(NaN, size(X0))
    n.ocp.XF        = XF
    n.ocp.XF_tol    = fill(NaN, size(XF))
    n.ocp.XL        = XL
    n.ocp.XU        = XU
    n.ocp.CL        = CL
    n.ocp.CU        = CU

    return n

end





function defineSolver!(n, kw)
    if typeof(kw)!=Dict{Symbol,Symbol}
      kw = Dict(kw)
    end

    # Default solver name is :Ipopt
    if haskey(kw, :name)
        n.s.ocp.solver.name = kw[:name]
    else
        @warn "No solver specified: using default solver (IPOPT)"
        n.s.ocp.solver.name = :Ipopt
    end

    if n.s.ocp.solver.name == :Ipopt
        n.ocp.mdl = Model(optimizer_with_attributes(Ipopt.Optimizer, _Ipopt_defaults...))
    else
        error("Solver $(n.s.ocp.solver.name) not defined")
    end
    
    return nothing

end













function OCPdef!(n::NLOpt{T}) where { T <: Number }

    # State variables
    varx = @variable(n.ocp.mdl, n.ocp.XL[i] <= x[j in 1:n.ocp.state.pts, i in 1:n.ocp.state.num] <= n.ocp.XU[i])
    n.r.ocp.x = varx

    # control variables
    varu = @variable(n.ocp.mdl,n.ocp.CL[i] <= u[j in 1:n.ocp.control.pts, i in 1:n.ocp.control.num] <= n.ocp.CU[i])
    n.r.ocp.u = varu

    # boundary constraints
    n.r.ocp.x0Con = []
    n.r.ocp.xfCon = []






    for st in 1:n.ocp.state.num
        if !isnan(n.ocp.X0[st]) # could have a bool for this
            n.r.ocp.x0Con = [ n.r.ocp.x0Con ; @constraint(n.ocp.mdl, x[1,st] == n.ocp.X0[st]) ]
                # fix(x[1, st], n.ocp.X0[st]; force = true)
        end
        if !isnan(n.ocp.XF[st])
            n.r.ocp.xfCon = [n.r.ocp.xfCon ; @constraint(n.ocp.mdl, x[end,st] == n.ocp.XF[st]) ]
        end
    end


    if n.s.ocp.integrationMethod == :ps
            n.r.ocp.dynCon = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]
            dynamics_expr = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]

            if n.s.ocp.finalTimeDV
                @variable(n.ocp.mdl, n.s.ocp.tfMin <= tf <=  n.s.ocp.tfMax)
                n.ocp.tf = tf
            end
            create_tV!(n)          # make a time vector

            for int in 1:n.ocp.Ni
                x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)

                # dynamics
                L = size(x_int)[1]-1;
                dx = Matrix{Any}(undef, L,n.ocp.state.num)
                for st in 1:n.ocp.state.num
                    dx[:,st] = DiffEq(n,x_int,u_int,L,st)
                end

                for st in 1:n.ocp.state.num # TODO consider multiplying X*D to reduce computations (i.e. remove this for loop for the states)
                    if n.s.ocp.integrationScheme == :lgrExplicit
                        dynamics_expr[int][:,st] = @expression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], sum(n.ocp.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.ocp.Nck[int]+1) - ((n.ocp.tf)/2)*dx[j,st]  )
                    elseif n.s.ocp.integrationScheme==:lgrImplicit
                        dynamics_expr[int][:,st] = @expression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], x_int[j+1,st] - x_int[1,st] - ((n.ocp.tf)/2)*sum(n.ocp.IMatrix[int][j,i]*dx[i,st] for i in 1:n.ocp.Nck[int]) )
                    end
                    for j in 1:n.ocp.Nck[int]
                        n.r.ocp.dynCon[int][j,st] = @constraint(n.ocp.mdl, 0. == dynamics_expr[int][j,st])
                    end
                end

                # additional constraints
                for num in 1:length(n.ocp.NLcon)
                    ch = addCon(n,x_int,u_int,L,num)
                    newConstraint!(n,ch,Symbol(string("ch",num))) # TODO could let the user name these
                end

            end

    elseif n.s.ocp.integrationMethod == :tm
        n.r.ocp.dynCon = Array{Any}(undef,n.ocp.N,n.ocp.state.num)

        if n.s.ocp.finalTimeDV
            @variable(n.ocp.mdl, n.s.ocp.tfMin <= tf <= n.s.ocp.tfMax)
            n.ocp.tf = tf
        end
        create_tV!(n) # make a time vector
        n.ocp.dt = n.ocp.tf / n.ocp.N * ones(n.ocp.N,)
        L = size(n.r.ocp.x)[1]
        #dx = [ DiffEq(n, n.r.ocp.x, n.r.ocp.u, L, st) for st in 1:n.ocp.state.num ]
        dx = Array{Any}(undef,L,n.ocp.state.num)
        for st in 1:n.ocp.state.num
          dx[:,st] = DiffEq(n,n.r.ocp.x,n.r.ocp.u,L,st)
        end
        if n.s.ocp.integrationScheme==:bkwEuler
          for st in 1:n.ocp.state.num
            n.r.ocp.dynCon[:,st] = @constraint(n.ocp.mdl, [j in 1:n.ocp.N], n.r.ocp.x[j+1,st] - n.r.ocp.x[j,st] ==  dx[j+1,st]*n.ocp.tf/(n.ocp.N) );
          end
        elseif n.s.ocp.integrationScheme==:trapezoidal
          for st in 1:n.ocp.state.num
            n.r.ocp.dynCon[:,st] = @constraint(n.ocp.mdl, [j in 1:n.ocp.N], n.r.ocp.x[j+1,st] - n.r.ocp.x[j,st] == 0.5*(dx[j,st] + dx[j+1,st])*n.ocp.tf/(n.ocp.N) )
          end
        else
            error("Not implemented yet")
        end

        # additional constraints
        for num in 1:length(n.ocp.NLcon)
            ch = addCon(n,n.r.ocp.x,n.r.ocp.u,L,num)
            newConstraint!(n,ch,Symbol(string("ch",num))) # TODO could let the user name these
        end
    end




    # save constraint data
    newConstraint!(n, n.r.ocp.x0Con , :x0_con )
    newConstraint!(n, n.r.ocp.xfCon , :xf_con )
    newConstraint!(n, n.r.ocp.dynCon, :dyn_con)

    # save the current working directory for navigation purposes
    n.r.mainDir = pwd()

    return nothing

end









function configure!(n::NLOpt{T}; kwargs... ) where { T <: Number }

    # Gather keyword arguments
    kw = Dict(kwargs)

    # Determine whether final time is a design variable (:finalTimeDV)
    # defaults to not a design variable
    if !haskey(kw,:finalTimeDV);n.s.ocp.finalTimeDV=false;
    else; n.s.ocp.finalTimeDV=get(kw,:finalTimeDV,0);
    end

    # Determine final time (:tf)
    if !haskey(kw,:tf) && !n.s.ocp.finalTimeDV
      error("\n If the final is not a design variable pass it as: (:tf=>Float64(some #)) \n
          If the final time is a design variable, indicate that as: (:finalTimeDV=>true)\n")
    elseif haskey(kw,:tf) && !n.s.ocp.finalTimeDV
      n.ocp.tf = Float64(get(kw,:tf,0))
    end



    # Integration Scheme (:integrationScheme)
    # Default = :lgrExplicit
    n.s.ocp.integrationScheme = get(kw, :integrationScheme, :bkwEuler)

    # Integration Method (:integrationMethod)
    if in( n.s.ocp.integrationScheme,  [ :lgrExplicit , :lgrImplicit ])

        # Use Pseudospectral Integration Method
        n.s.ocp.integrationMethod = :ps

        if haskey(kw, :N)
            # error(":N is not an appropriate keyword argument for :ps method :$(n.s.ocp.integrationScheme), use :Nck")
        else

            n.ocp.Nck = get(kw, :Nck, [10 , 10 , 10 , 10] )

            if any(n.ocp.Nck .< 0)
                error(":Nck must be > 0")
            end

            # Determine number of points
            n.ocp.Ni            = length(n.ocp.Nck)
            n.ocp.state.pts     = sum(n.ocp.Nck) + 1
            n.ocp.control.pts   = sum(n.ocp.Nck)
            n.ocp.Nck_full      = [ 0 ; cumsum(n.ocp.Nck .+ 1) ]
            n.ocp.Nck_cum       = [ 0 ; cumsum(n.ocp.Nck) ]

            # Initialize node data
            taus_and_weights    = [ gaussradau(n.ocp.Nck[int]) for int in 1:n.ocp.Ni ]
            n.ocp.tau           = [taus_and_weights[int][1] for int in 1:n.ocp.Ni]
            n.ocp.w             = [taus_and_weights[int][2] for int in 1:n.ocp.Ni]

            # Create Intervals
            createIntervals!(n)

            #? Create DMatrix
            DMatrix!(n)

        end


    elseif in( n.s.ocp.integrationScheme,  [ :trapezoidal, :bkwEuler ] )

        # Use Trapezoidal Method
        n.s.ocp.integrationMethod = :tm

        if haskey(kw, :Nck)
            error(":Nck is not an appropriate keyword argument for :tm method :$(n.s.ocp.integrationScheme), use :N")
        else

            # Default number of points is 100
            n.ocp.N = get(kw, :N, 100)

            # Set number of integration points for state and control
            n.ocp.state.pts = n.ocp.N + 1
            n.ocp.control.pts = n.ocp.N + 1

        end

    else
        error("The :integrationScheme that you specified ($(n.s.ocp.integrationScheme)) is not currently implemeted.")
    end
    n.ocp.mXL = falses(n.ocp.state.num)
    n.ocp.mXU = falses(n.ocp.state.num)
    #n.ocp.XL_var = Vector{Vector{JuMP.JuMPTypes}}([ [ @variable(JuMP.Model(), tmp) for i = 1:n.ocp.state.num ] for j = 1:n.ocp.state.pts] )
    #n.ocp.XU_var = Vector{Vector{JuMP.JuMPTypes}}([ [ @variable(JuMP.Model(), tmp) for i = 1:n.ocp.state.num ] for j = 1:n.ocp.state.pts] )
    n.ocp.XL_var = Matrix{Float64}(undef,n.ocp.state.num,n.ocp.state.pts)
    n.ocp.XU_var = Matrix{Float64}(undef,n.ocp.state.num,n.ocp.state.pts)

    # Define Solver Settings
    #TODO Change here to accomodate MadNLP
    defineSolver!(n, Dict(get(kw, :solverSettings, (:name => :Ipopt) )))

    # Definie Optimal Control Problem
    OCPdef!(n)

    return nothing

end

