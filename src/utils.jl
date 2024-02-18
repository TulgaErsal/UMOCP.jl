include("types.jl")
include("math.jl")
include("setup.jl")
include("ps.jl")
include("diffeq.jl")

function intervals(n::NLOpt,int,x,u)

  if typeof(x[1,1]) == JuMP.VariableRef

    # States
    x_int = Matrix{JuMP.VariableRef}(undef, length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.state.num)
    for st in 1:n.ocp.state.num # +1 adds the DV in the next interval
      x_int[:,st] = x[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,st]
    end

    # Controls
    if int!=n.ocp.Ni
      u_int = Matrix{JuMP.VariableRef}(undef, length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.control.num)
    else                    # -1 -> removing control in last mesh interval
      u_int = Matrix{JuMP.VariableRef}(undef, length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]-1),n.ocp.control.num)
    end
    for ctr in 1:n.ocp.control.num
      if int!=n.ocp.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,ctr]
      else
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1],ctr]
      end
    end

  else
    # states
    x_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.state.num);
    for st in 1:n.ocp.state.num # +1 adds the DV in the next interval
      x_int[:,st] = x[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,st]
    end

    # controls
    if int!=n.ocp.Ni
      u_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.control.num)
    else                    # -1 -> removing control in last mesh interval
      u_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]-1),n.ocp.control.num)
    end
    for ctr in 1:n.ocp.control.num
      if int!=n.ocp.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,ctr]
      else
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1],ctr]
      end
    end
  end

  return x_int,u_int
end


function initState(numStates)::State
  s = State()
  s.num = numStates
  s.name = [Symbol("x$i") for i in 1:numStates]
  s.description = [String("x$i") for i in 1:numStates]
  return s
end
function initControl(numControls)::Control
  c = Control()
  c.num = numControls
  c.name = [Symbol("u$i") for i in 1:numControls]
  c.description = [String("u$i") for i in 1:numControls]
  return c
end

function states!(n::NLOpt,names;descriptions=[])

  if length(names)!=n.ocp.state.num
    error("\n Check size of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.state.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.state.num
    if names[i]==:xxx
      error("xxx is OFF limits for a state name; please choose something else. \n")
    end
  end

   n.ocp.state.name = names
   if !isempty(descriptions)
     n.ocp.state.description = descriptions
   end
  return nothing
end


function controls!(n::NLOpt,names;descriptions=[])
  if length(names)!=n.ocp.control.num
    error("\n Check sizes of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.control.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.control.num
    if names[i]==:uuu
      error("uuu is OFF limits for a control name; please choose something else. \n")
    end
   end

  n.ocp.control.name = names
  if !isempty(descriptions)
    n.ocp.control.description = descriptions
  end
  return nothing
end


function defineTolerances!(n::NLOpt;
                          X0_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,),
                          XF_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,))
  # TODO error checking, if the user does not pass tolerances etc.
  n.ocp.X0_tol = X0_tol
  n.ocp.XF_tol = XF_tol
  return nothing
end


function create_tV!(n::NLOpt)

  if n.s.ocp.integrationMethod==:ps
    # create mesh points, interval size = tf_var/Ni
    tm = @expression(n.ocp.mdl, [idx=1:n.ocp.Ni+1], (idx-1)*n.ocp.tf/n.ocp.Ni)
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts = [Array{Any}(undef, n.ocp.Nck[int]+1,) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      ts[int][1:end-1] = @expression(n.ocp.mdl,[j=1:n.ocp.Nck[int]], (tm[int+1]-tm[int])/2*n.ocp.tau[int][j] +  (tm[int+1]+tm[int])/2);
      ts[int][end] = @expression(n.ocp.mdl, n.ocp.tf/n.ocp.Ni*int) # append +1 at end of each interval
    end
    tt1 = [idx for tempM in ts for idx = tempM[1:end-1]]
    tmp = [tt1;ts[end][end]]
    n.ocp.tV = @expression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
  else
    # create vector with the design variable in it
    t = Array{Any}(undef, n.ocp.N+1,1)
    tm = @expression(n.ocp.mdl, [idx=1:n.ocp.N], n.ocp.tf/n.ocp.N*idx)
    tmp = [0;tm]
    n.ocp.tV = @expression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
  end
  return nothing
end


function initConstraint!(n::NLOpt)
    if n.r.ocp.constraint === nothing
        n.r.ocp.constraint = Constraint()
    end
    return nothing
end

"""
"""
function newConstraint!(n::NLOpt, handle, name::Symbol)
    initConstraint!(n)
    push!(n.r.ocp.constraint.name,name)
    push!(n.r.ocp.constraint.handle,handle)

    return nothing
end


function integrate!(n::NLOpt,V::Expr)
  if n.s.ocp.integrationMethod==:ps
    integral_expr = [Array{Any}(undef, n.ocp.Nck[int]) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)
      L = size(x_int)[1]-1
      integral_expr[int][:] = NLExpr(n,V,x_int,u_int,L)
    end
    @expression(n.ocp.mdl, temp[int=1:n.ocp.Ni], (n.ocp.tf-n.ocp.t0)/2*sum(n.ocp.ws[int][j]*integral_expr[int][j] for j = 1:n.ocp.Nck[int]) )
    expression = @expression(n.ocp.mdl, sum(temp[int] for int = 1:n.ocp.Ni))
  elseif n.s.ocp.integrationMethod==:tm
    L = size(n.r.ocp.x)[1]
    temp = NLExpr(n,V,n.r.ocp.x,n.r.ocp.u,L);
    if n.s.ocp.integrationScheme==:bkwEuler
      # NOTE integration this way does not penalize the first control
      expression = @expression(n.ocp.mdl, sum(temp[j+1]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
      #expression = @NLexpression(n.ocp.mdl, sum(temp[j]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    elseif n.s.ocp.integrationScheme ==:trapezoidal || n.s.ocp.integrationScheme == :mpcol
      expression = @expression(n.ocp.mdl, sum(0.5*(temp[j]+temp[j+1])*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    else
      error("\n $(n.s.ocp.integrationScheme) not defined in integrationSchemes\n")
    end
  else
    error("\n $(n.s.ocp.integrationMethod) not defined in integrationMethods \n")
  end
  return expression
end



function OCPoptimize!(n::NLOpt)

    status = JuMP.optimize!(n.ocp.mdl)

    if !n.s.ocp.cacheOnly
        n.r.ocp.Terminalstatus = termination_status(n.ocp.mdl)
        n.r.ocp.status = RetrieveSolveStatus(n)
        n.r.ocp.tSolve = solve_time(n.ocp.mdl)
        n.r.ocp.objVal = objective_value(n.ocp.mdl)
        n.r.ocp.iterNum = barrier_iterations(n.ocp.mdl)    # possible iteration number for a higher level algorithm
        n.r.ocp.evalNum = n.r.ocp.evalNum + 1
        postProcess!(n)      # temporarily save data
    end
    WarmStart(n)
    return nothing
end


function postProcess!(n::NLOpt; kwargs...)

    kw = Dict(kwargs)


    if n.s.ocp.save
        opt2dfs!(n)
    end

    if ((n.r.ocp.status==:Optimal) || (n.r.ocp.status==:UserLimit))
        if n.s.ocp.integrationMethod == :ps
            if n.s.ocp.finalTimeDV
                t = [scale_tau(n.ocp.ts[int],0.0,value.(n.ocp.tf)) for int in 1:n.ocp.Ni]     # scale time from [-1,1] to [t0,tf]
            else
                t = [scale_tau(n.ocp.ts[int],0.0,n.ocp.tf) for int in 1:n.ocp.Ni]
            end
            n.r.ocp.tctr = [idx for tempM in t for idx = tempM[1:end-1]] .+ n.ocp.t0
            n.r.ocp.tst = [n.r.ocp.tctr; t[end][end] .+ n.ocp.t0]
            # TODO check the above line... is t0 getting added on twice?
        elseif n.s.ocp.integrationMethod == :tm
            if n.s.ocp.finalTimeDV
                n.r.ocp.tctr = append!([0.0],cumsum(value.(n.ocp.dt))) .+ n.ocp.t0
            else
                n.r.ocp.tctr = append!([0.0],cumsum(n.ocp.dt)) .+ n.ocp.t0
            end
            n.r.ocp.tst = n.r.ocp.tctr
        end

        stateDataExists = false
        if (n.r.ocp.status==:Optimal) || (n.r.ocp.status==:UserLimit)
            stateDataExists = true
            n.r.ocp.X = zeros(Float64,n.ocp.state.pts,n.ocp.state.num)
            n.r.ocp.U = zeros(Float64,n.ocp.control.pts,n.ocp.control.num)
            for st in 1:n.ocp.state.num
                n.r.ocp.X[:,st] = value.(n.r.ocp.x[:,st])
            end
            for ctr in 1:n.ocp.control.num
                n.r.ocp.U[:,ctr] = value.(n.r.ocp.u[:,ctr])
            end
        else
            @warn "The solution is not Optimal \n"
        end

        if n.s.ocp.save && stateDataExists
            if n.s.ocp.InternalLogging
                push!(n.r.ocp.dfs,dvs2dfs(n))
            else
                n.r.ocp.dfs = Vector{DataFrame}()
                push!(n.r.ocp.dfs,dvs2dfs(n))
            end
        end
    end

    return nothing
end


function opt2dfs!(n::NLOpt;kwargs...)

    kw = Dict(kwargs)

    if !haskey(kw,:statusUpdate)
        statusUpdate = false
    else
        statusUpdate = get(kw,:statusUpdate,0)
    end

    # make sure that the feildnames are initialized
    if isempty(n.r.ocp.dfsOpt)
        n.r.ocp.dfsOpt = DataFrame(tSolve = [], objVal = [], status = [], iterNum = [], evalNum = [])
    end

    if !statusUpdate
        push!(n.r.ocp.dfsOpt[!, :tSolve], n.r.ocp.tSolve)
        push!(n.r.ocp.dfsOpt[!, :objVal], n.r.ocp.objVal)
        push!(n.r.ocp.dfsOpt[!, :status], n.r.ocp.status)
    else  # TODO consider removing and cleaning this up
        push!(n.r.ocp.dfsOpt[!, :tSolve], NaN)
        push!(n.r.ocp.dfsOpt[!, :objVal], NaN)
        if statusUpdate && (typeof(n.r.ocp.status)==Symbol)
            push!(n.r.ocp.dfsOpt[!, :status], n.r.ocp.status)
        else
            push!(n.r.ocp.dfsOpt[!, :status], NaN)
        end
    end
    push!(n.r.ocp.dfsOpt[!, :evalNum], n.r.ocp.evalNum-1)

    return nothing
end


function dvs2dfs(n::NLOpt)

    dfs = DataFrame()
    dfs[!, :t] = n.r.ocp.tst
    for st in 1:n.ocp.state.num
        dfs[!, n.ocp.state.name[st]] = n.r.ocp.X[:,st]
    end
    for ctr in 1:n.ocp.control.num
        if n.s.ocp.integrationMethod==:tm
            dfs[!, n.ocp.control.name[ctr]] = n.r.ocp.U[:,ctr]
        else
            dfs[!, n.ocp.control.name[ctr]] = [n.r.ocp.U[:,ctr];NaN]
        end
    end

    return dfs
end


function WarmStart(n::NLOpt)
  flag = false
  try get_attribute(n.ocp.mdl, "warm_start_init_point")
    if get_attribute(n.ocp.mdl, "warm_start_init_point") == "yes"
      flag = true
    else
      flag = false
    end
  catch
    flag = false
  end
  if flag == true
    set_start_value.(n.r.ocp.x, n.r.ocp.X)
    set_start_value.(n.r.ocp.u, n.r.ocp.U)
  end

end


function RetrieveSolveStatus(n::NLOpt)
  status = termination_status(n.ocp.mdl)
   
  SolvingStatus = [:Optimal, :UserLimit, :InFeasible]
  OptimalList = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
  LimitList = [MOI.ITERATION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.SOLUTION_LIMIT, MOI.MEMORY_LIMIT, MOI.OBJECTIVE_LIMIT, MOI.NORM_LIMIT, MOI.OTHER_LIMIT ]
  
  if status ∈ OptimalList
      return SolvingStatus[1]
  elseif status ∈ LimitList
      return SolvingStatus[2]
  else
      return SolvingStatus[3]
  end
end