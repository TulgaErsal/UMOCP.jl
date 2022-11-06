function try_import(name::Symbol)
    try
        @eval using $name
        return true
    catch e
    end
end


function intervals(n::NLOpt,int,x,u)

  if typeof(x[1,1]) == JuMP.VariableRef

    # States
    x_int = Array{JuMP.VariableRef}(length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.state.num)
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
  if !n.f.ocp.defined
    error("\n call define() before calling states!() \n")
  end
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
  if !n.f.ocp.defined
    error("\n call define() before calling controls!() \n")
  end
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
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.Ni+1], (idx-1)*n.ocp.tf/n.ocp.Ni)
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts = [Array{Any}(undef, n.ocp.Nck[int]+1,) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      ts[int][1:end-1] = @NLexpression(n.ocp.mdl,[j=1:n.ocp.Nck[int]], (tm[int+1]-tm[int])/2*n.ocp.tau[int][j] +  (tm[int+1]+tm[int])/2);
      ts[int][end] = @NLexpression(n.ocp.mdl, n.ocp.tf/n.ocp.Ni*int) # append +1 at end of each interval
    end
    tt1 = [idx for tempM in ts for idx = tempM[1:end-1]]
    tmp = [tt1;ts[end][end]]
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
  else
    # create vector with the design variable in it
    t = Array{Any}(undef, n.ocp.N+1,1)
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.N], n.ocp.tf/n.ocp.N*idx)
    tmp = [0;tm]
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
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
    # TODO test if this is the test that is making sure constraints can be entered in results in postProcess
    #@info "newConstraint!: testing vcat:                        $(vcat(n.r.ocp.constraint.handle, handle))"

    push!(n.r.ocp.constraint.name,name)
    #n.r.ocp.constraint.handle = vcat(n.r.ocp.constraint.handle, handle)
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
    @NLexpression(n.ocp.mdl, temp[int=1:n.ocp.Ni], (n.ocp.tf-n.ocp.t0)/2*sum(n.ocp.ws[int][j]*integral_expr[int][j] for j = 1:n.ocp.Nck[int]) )
    expression = @NLexpression(n.ocp.mdl, sum(temp[int] for int = 1:n.ocp.Ni))
  elseif n.s.ocp.integrationMethod==:tm
    L = size(n.r.ocp.x)[1]
    temp = NLExpr(n,V,n.r.ocp.x,n.r.ocp.u,L);
    if n.s.ocp.integrationScheme==:bkwEuler
      # NOTE integration this way does not penalize the first control
      expression = @NLexpression(n.ocp.mdl, sum(temp[j+1]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
      #expression = @NLexpression(n.ocp.mdl, sum(temp[j]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    elseif n.s.ocp.integrationScheme ==:trapezoidal || n.s.ocp.integrationScheme == :mpcol
      expression = @NLexpression(n.ocp.mdl, sum(0.5*(temp[j]+temp[j+1])*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    else
      error("\n $(n.s.ocp.integrationScheme) not defined in integrationSchemes\n")
    end
  else
    error("\n $(n.s.ocp.integrationMethod) not defined in integrationMethods \n")
  end
  return expression
end



function NLoptimize!(n::NLOpt; Iter::Int=0)

    status = JuMP.optimize!(n.ocp.mdl)

    if !n.s.ocp.cacheOnly
        n.r.ocp.status = termination_status(n.ocp.mdl)
        n.r.ocp.tSolve = solve_time(n.ocp.mdl)
        n.r.ocp.objVal = objective_value(n.ocp.mdl)
        n.r.ocp.iterNum = Iter    # possible iteration number for a higher level algorithm
        n.r.ocp.evalNum = n.r.ocp.evalNum + 1
        postProcess!(n)      # temporarily save data
    end

    return nothing
end

function evalConstraints!(n::NLOpt)

    n.r.ocp.constraint.value = []   # reset values
    n.r.ocp.constraint.nums = []
    s=1

    for i = 1:length(n.r.ocp.constraint.handle)
        if n.r.ocp.constraint.name[i] == :dyn_con  # state constraits
            dfs=Vector{DataFrame}(undef, n.ocp.state.num)
            con=DataFrame(step=1)
            l=0
            for st in 1:n.ocp.state.num
                if n.s.ocp.integrationMethod==:ps
                    temp= [ dual.(n.r.ocp.constraint.handle[i][int][:,st]) for int in 1:n.ocp.Ni ]
                    vals = [ idx for tempM in temp for idx in tempM ]
                    dfs[st] = DataFrame(step=1:sum(n.ocp.Nck) ; Dict(n.ocp.state.name[st] => vals)...)
                    l += length(vals)
                else
                    dfs[st] = DataFrame(step=1:n.ocp.N; Dict(n.ocp.state.name[st] => dual.(n.r.ocp.constraint.handle[i][:,st]))...)
                    l += length(n.r.ocp.constraint.handle[i][:,st])
                end

                con = ( st == 1 ? dfs[st] : innerjoin(con, dfs[st], on=:step) )

            end
        else
            S=0
            try
                S=JuMP.size(n.r.ocp.constraint.handle[i])
            catch
                error("\n For now, the constraints cannot be in this form: \n
                con=@NLconstraint(mdl,n.r.ocp.u[1,1]==param); \n
                Write it in array form: \n
                con=@NLconstraint(mdl,[i=1],n.r.ocp.u[i,1]==param); \n")
            end
            if length(S) == 1
                if typeof(n.r.ocp.constraint.handle[i]) == Vector{NonlinearConstraintRef{ScalarShape}} || typeof(n.r.ocp.constraint.handle[i]) == Vector{Any}
                    con = DataFrame(step=1:length(n.r.ocp.constraint.handle[i]); Dict(n.r.ocp.constraint.name[i] => dual.(n.r.ocp.constraint.handle[i][:]))...)
                elseif typeof(n.r.ocp.constraint.handle[i]) == JuMP.Containers.DenseAxisArray{NonlinearConstraintRef{ScalarShape}, 1, Tuple{UnitRange{Int64}}, Tuple{JuMP.Containers._AxisLookup{Tuple{Int64, Int64}}}}
                    handle_key = axes(n.r.ocp.constraint.handle[i])[1]
                    dual_con = Vector{NonlinearConstraintRef{ScalarShape}}[]
                    for key in handle_key
                        dual_con = [dual_con; n.r.ocp.constraint.handle[6][key]]
                    end
                    con = DataFrame(step=1:length(n.r.ocp.constraint.handle[i]); Dict(n.r.ocp.constraint.name[i] => dual.(dual_con))...)
                  else
                    error("\n This is not the correct of Constraint!")
                end
                l=S[1]


            elseif length(S) == 2
                dfs=Vector{DataFrame}(S[1])
                con=DataFrame(step=1)
                for idx in 1:S[1]
                    try
                        dfs[idx] = DataFrame(step=1:S[2]; Dict(n.r.ocp.constraint.name[i] => dual.(n.r.ocp.constraint.handle[i][idx,:]))...)
                    catch
                        dfs[idx] = DataFrame(step=1:S[2]; Dict(n.r.ocp.constraint.name[i] => NaN)...) # fix for case where all of the states are not being constrainted, but some are within some XF_tol
                    end
                    if idx==1;
                        con = dfs[idx]
                    else
                        con = join(con,dfs[idx],on=:step,makeunique=true)
                    end
                end

                l = S[1] * S[2]
            end
        end
        f = s + l - 1
        num = (i, n.r.ocp.constraint.name[i], "length = $l", string("indices in g(x) = "), (s, f))
        push!(n.r.ocp.constraint.nums,num)
        push!(n.r.ocp.constraint.value,con)
        s = f + 1
    end
    return nothing
end

function postProcess!(n::NLOpt; kwargs...)

    kw = Dict(kwargs)
    # check to see if the user is initializing while compensating for control delay
    if !haskey(kw,:Init)
        Init = false
    else
        Init = get(kw,:Init,0)
    end

    if n.s.ocp.save
        opt2dfs!(n)
    end

    # even if n.r.ocp.status==:Infeasible try to get solution. For the case that user may want to look at results to see where constraints where violated
    # in this case set =>  n.s.ocp.evalConstraints = true
    # http://jump.readthedocs.io/en/latest/refmodel.html#solve-status
    if !Init && (n.s.ocp.evalConstraints || ((n.r.ocp.status==MOI.OPTIMAL) || (n.r.ocp.status==MOI.LOCALLY_SOLVED)))
        if n.s.ocp.integrationMethod == :ps
            if n.s.ocp.finalTimeDV
                t = [scale_tau(n.ocp.ts[int],0.0,value.(n.ocp.tf)) for int in 1:n.ocp.Ni]     # scale time from [-1,1] to [t0,tf]
            else
                t = [scale_tau(n.ocp.ts[int],0.0,n.ocp.tf) for int in 1:n.ocp.Ni]
            end
            n.r.ocp.tctr = [idx for tempM in t for idx = tempM[1:end-1]] .+ value(n.ocp.t0)
            n.r.ocp.tst = [n.r.ocp.tctr; t[end][end] .+ value(n.ocp.t0)]
            # TODO check the above line... is t0 getting added on twice?
        elseif n.s.ocp.integrationMethod == :tm
            if n.s.ocp.finalTimeDV
                n.r.ocp.tctr = append!([0.0],cumsum(value.(n.ocp.dt))) .+ value(n.ocp.t0)
            else
                n.r.ocp.tctr = append!([0.0],cumsum(n.ocp.dt)) .+ value(n.ocp.t0)
            end
            n.r.ocp.tst = n.r.ocp.tctr
        end

        stateDataExists = false
        if n.r.ocp.status == MOI.OPTIMAL || n.r.ocp.status == MOI.LOCALLY_SOLVED
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
        if n.s.ocp.evalConstraints && n.r.ocp.status != INFEASIBLE && n.r.ocp.status != LOCALLY_INFEASIBLE && n.r.ocp.status != OTHER_ERROR # note may want to remove the && arg
            evalConstraints!(n)
            # # TODO: Evaluated constraints in n.r.ocp.constraints
            # if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps
            #     L1 = 0       # find index where dynamics constraints start
            #     for i in 1:length(n.r.ocp.constraint.name)
            #         if n.r.ocp.constraint.name[i] == :dyn_con
            #             L1 = n.r.ocp.constraint.nums[i][end][1]
            #         end
            #     end
            #     mpb = JuMP.internalmodel(n.ocp.mdl)
            #     c = MathProgBase.getconstrduals(mpb)
            #     # NOTE for now since costates are not defined for :tm methods, n.r.ocp.CS is in a different format than n.r.ocp.X
            #     # in the future if costates are defined for :tm methods this can be changed
            #     n.r.ocp.CS = [[zeros(Float64,n.ocp.Nck[int]) for int in 1:n.ocp.Ni] for st in 1:n.ocp.state.num]
            #     for int in 1:n.ocp.Ni
            #         b = 0
            #         for st in 1:n.ocp.state.num
            #             a = L1 + n.ocp.Nck[int]*(st-1)  # n.ocp.Nck[int]*(st-1) adds on indices for each additional state within the interval
            #             b = a + n.ocp.Nck[int] - 1      # length of the state within this interval
            #             n.r.ocp.CS[st][int] = -c[a:b]./n.ocp.ws[int]
            #         end
            #         L1 = b + 1 # adds indicies due to a change in the interval
            #     end
            # end
        end

        if n.s.ocp.save && stateDataExists
            push!(n.r.ocp.dfs,dvs2dfs(n))
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

    # if n.s.ocp.integrationMethod == :ps && n.s.ocp.evalConstraints
    #     CS_vector = Matrix{Float64}(undef, n.ocp.state.pts, n.ocp.state.num)
    #     for st in 1:n.ocp.state.num # states
    #         temp = [n.r.ocp.CS[st][int][1:end,:] for int in 1:n.ocp.Ni]
    #         CS_vector[1:end-1,st] = [idx for tempM in temp for idx=tempM]
    #         CS_vector[end,st] = NaN
    #     end
    #     for st in 1:n.ocp.state.num # states
    #         dfs[!, Symbol(n.ocp.state.name[st],:cs)] = CS_vector[:,st]
    #     end
    # end

    return dfs
end
