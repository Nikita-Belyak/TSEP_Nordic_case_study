function single_level_problem_generation(ip::initial_parameters)    
    # Defining single-level model
    single_level_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    #single_level_problem = Model(dual_optimizer(() -> Gurobi.Optimizer(GRB_ENV)))
    #single_level_problem = Model(() -> Ipopt.Optimizer())
    #set_optimizer_attribute(single_level_problem, "OutputFlag", 0)
    set_optimizer_attribute(single_level_problem, "NonConvex", 2)
    #set_optimizer_attribute(single_level_problem, "InfUnbdInfo", 1)
    #set_optimizer_attribute(single_level_problem, "Presolve", 0)
    #set_optimizer_attribute(single_level_problem, "IntFeasTol", 1E-9)
    #set_optimizer_attribute(single_level_problem, "FeasibilityTol", 1E-9)
    #set_optimizer_attribute(single_level_problem, "FeasibilityTol", 1E-9)
    #set_optimizer_attribute(single_level_problem, "NumericFocus", 2)

    ## VARIABLES

    # Conventional generation related variable
    @variable(single_level_problem, g_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # VRES generation related variable
    @variable(single_level_problem, g_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)

    # hydro power generation related variable
    @variable(single_level_problem, g_hydro[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod] >= 0)

    # Consumption realted variable
    @variable(single_level_problem, q[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes] >= 0)

    # Energy transmission realted variable
    @variable(single_level_problem, f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes]) # >= 0)

    # Transmission capacity expansion realted variable
    @variable(single_level_problem, l_plus[ 1:ip.num_nodes, 1:ip.num_nodes]>=0)

    #addition 
    @variable(single_level_problem, 1 <= y[ 1:ip.num_nodes, 1:ip.num_nodes] <= 4, Int)

    # Conventional energy capacity expansion realted variable
    @variable(single_level_problem, g_conv_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv]>=0)

    # VRES capacity expansion realted variable
    @variable(single_level_problem, g_VRES_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES]>=0)

    ## OBJECTIVE
    @objective(single_level_problem, Max, 
            sum(
                sum( ip.scen_prob[s] * 

                    (ip.id_intercept[s,t,n]*q[s,t,n] 
                    - 0.5* ip.id_slope[s,t,n]/ip.time_periods[t]* q[s,t,n]^2
                    
                    - sum( 
                        (ip.conv.operational_costs[n,i,e] 
                        + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                    for i in 1:ip.num_prod, e in 1:ip.num_conv)

                    #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                    #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                    #- 1E-4 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                    )

                    

                for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
                
                -sum( 

                    sum( 
                        ip.vres.maintenance_costs[n,i,r]*
                        (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                        + 
                        ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                    for r in 1:ip.num_VRES)
                    
                    +
                        
                    sum( 
                        ip.conv.maintenance_costs[n,i,e]*
                        (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                        + 
                        ip.conv.investment_costs[n,i,e]* g_conv_plus[n,i,e]
                    for e in 1:ip.num_conv) 
                         
                for i in 1:ip.num_prod)
                
                - sum(
                    ip.transm.maintenance_costs[n,m]*
                    (ip.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m])
                    +
                    ip.transm.investment_costs[n,m] * 0.5 * l_plus[n,m]
                for m in 1:ip.num_nodes)
            
            for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )
    
    ## CONSTRAINTS

    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        q[s,t,n]/scaling_factor  - sum( 
                    sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
                    + 
                    sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES)
                    +
                    sum(g_hydro[s,t,n,i]) 
                    for i = 1:ip.num_prod) / scaling_factor 
                + sum( f[s,t,n,m] for m in n+1:ip.num_nodes)/scaling_factor - sum(f[s,t,m,n] for m in 1:n-1)/scaling_factor #- sum( 0.98*f[s,t,m,n] for m in 1:n-1)
        == 0
    )

    #addition 
    @constraint(single_level_problem, [n in 1:ip.num_nodes, m in 1:ip.num_nodes ], l_plus[n,m] == 3000*(y[n,m]-1)) #

    
    # Transmission bounds 
    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor <= 0
    )

    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor + ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor >= 0
    )
    
    # Primal feasibility for the transmission 
    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        f[s,t,n,m]/scaling_factor == 0 
    )

    #@constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        #f[s,t,n,m] - sum( 
            #sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
           #  + 
           # sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
           # for i = 1:ip.num_prod) - sum(f[s,t,m1,n] for m1 in 1:ip.num_nodes)<= 0
    #)

    # Generation budget limits
    @constraint(single_level_problem, [i = 1:ip.num_prod],
        sum(ip.vres.investment_costs[n,i,r] * g_VRES_plus[n,i,r] for  r in 1:ip.num_VRES, n in 1:ip.num_nodes)/scaling_factor 
        + sum(ip.conv.investment_costs[n,i,e] * g_conv_plus[n,i,e] for  e in 1:ip.num_conv, n in 1:ip.num_nodes)/scaling_factor <= ip.gen_budget[i]/scaling_factor
    )

    # Conventional generation bounds
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor - ip.time_periods[t]*(ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )

    # VRES generation bounds 
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        g_VRES[s,t,n,i,r]/scaling_factor - ip.time_periods[t]* ip.vres.availability_factor[s,t,n,r]*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r])/scaling_factor <= 0
    )

    # Hydro power generation bounds 
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod],
        g_hydro[s,t,n,i]/scaling_factor - (ip.hydro[s,t,n,i])/scaling_factor <= 0
    )

    # primal feasibility for the transmission expansion investments
    #@constraint(single_level_problem, [n in 1:ip.num_nodes],
        #l_plus[n,n] == 0
    #)

    @constraint(single_level_problem, [n in 1:ip.num_nodes, m in 1:ip.num_nodes],
        l_plus[n,m]/scaling_factor - l_plus[m,n]/scaling_factor == 0
    )

    # Primal feasibility for the transmission expansion 2 (budget related)
    @constraint(single_level_problem,
        sum(ip.transm.investment_costs[n,m] * 0.5 * l_plus[n,m] for n in 1:ip.num_nodes, m in 1:ip.num_nodes)/scaling_factor - ip.transm.budget_limit/scaling_factor <= 0 
    )



    # Maximum ramp-up rate for conventional units
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
       g_conv[s,t,n,i,e]/scaling_factor - (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor - ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )

    # Maximum ramp-down rate for conventional units
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor - g_conv[s,t,n,i,e]/scaling_factor - ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )
    
    return single_level_problem
end

function bi_level_problem_generation_new(ip::initial_parameters, market::String)    
    # Defining single-level model
    bi_level_problem = BilevelModel(() -> Gurobi.Optimizer(), mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = 1E+10, dual_big_M = 1E+10))
    #bi_level_problem = Model(() -> Ipopt.Optimizer())
    #set_optimizer_attribute(bi_level_problem, "OutputFlag", 0)
    #set_optimizer_attribute(bi_level_problem, "NonConvex", 2)
    #set_optimizer_attribute(single_level_problem, "InfUnbdInfo", 1)
    #set_optimizer_attribute(bi_level_problem, "Presolve", 0)
    #set_optimizer_attribute(bi_level_problem, "IntFeasTol", 1E-9)
    

    set_optimizer_attribute(bi_level_problem, "LogFile", src_link * "/2022 data/optimisation_results/"*Dates.format(Dates.now(),"dd_u_yyyy_HH_MM_SS")*".log" )
    set_optimizer_attribute(bi_level_problem, "Threads", Sys.CPU_THREADS)
    #set_optimizer_attribute(bi_level_problem, "MIPFocus", 3)


    #set_optimizer_attribute(bi_level_problem, "IntFeasTol", 1E-9)
    #set_optimizer_attribute(bi_level_problem, "IntFeasTol", 1E-9)
    
    #set_optimizer_attribute(bi_level_problem, "FeasibilityTol", 1E-9)
    #set_optimizer_attribute(bi_level_problem, "NumericFocus", 3)
    #set_optimizer_attribute(bi_level_problem, "BarCorrectors", 2000000000)
    #set_optimizer_attribute(bi_level_problem, "BarQCPConvTol", 0.0)
    #set_optimizer_attribute(bi_level_problem, "BarConvTol", 0.0)
    #set_optimizer_attribute(bi_level_problem, "BarHomogeneous", 1)
    #set_optimizer_attribute(bi_level_problem, "Method", 4)
    #set_optimizer_attribute(bi_level_problem, "MIPFocus", 3)
 

    #CPLEX

    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Strategy_Probe", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_Cliques", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_Covers", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_Disjunctive", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_LiftProj", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_LocalImplied", 3)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_BQP", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_FlowCovers", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_PathCut", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_Gomory", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_GUBCovers", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_MIRCut", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cut_MCFCut", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_RLT", 2)
    #set_optimizer_attribute(bi_level_problem, "CPXPARAM_MIP_Cuts_ZeroHalfCut", 2)

    #set_optimizer_attribute(bi_level_problem, "CPX_PARAM_VARSEL", 4)

    #set_optimizer_attribute(bi_level_problem, "CPX_PARAM_VARSEL", 3)

    ## UPPER LEVEL VARIABLES

    # Transmission capacity expansion realted variable
    @variable(Upper(bi_level_problem), l_plus[ 1:ip.num_nodes, 1:ip.num_nodes]>=0)

    ## LOWER LEVEL VARIABLES 

    # Conventional generation related variable
    @variable(Lower(bi_level_problem), g_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # VRES generation related variable
    @variable(Lower(bi_level_problem), g_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)
    
    # hydro power generation related variable
    @variable(Lower(bi_level_problem), g_hydro[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod] >= 0)

    # Consumption realted variable
    @variable(Lower(bi_level_problem), q[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes]>= 0)

    # Energy transmission realted variable
    @variable(Lower(bi_level_problem), f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes] ) #>= 0)

    # Conventional energy capacity expansion realted variable
    @variable(Lower(bi_level_problem), g_conv_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv]>=0)

    # VRES capacity expansion realted variable
    @variable(Lower(bi_level_problem), g_VRES_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES]>=0)  

    ## UPPER LEVEL OBJECTIVE

    @objective(Upper(bi_level_problem), Max, 
        sum(
            sum( ip.scen_prob[s] * 

                (ip.id_intercept[s,t,n]*q[s,t,n] 
                - 0.5* ip.id_slope[s,t,n]/ip.time_periods[t]* q[s,t,n]^2

                #- (market == "perfect" ? 0 : 
                    #(0.5* ip.id_slope[s,t,n] * sum( 
                    #    (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2
                    #for i in 1:ip.num_prod))
                    #)
            
                - sum( 
                    (ip.conv.operational_costs[n,i,e] 
                    + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                for i in 1:ip.num_prod, e in 1:ip.num_conv)

                #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                #- 1E-1 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                )

            for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
        
            -sum( 

                sum( 
                    ip.vres.maintenance_costs[n,i,r]*
                    (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                    + 
                    ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                for r in 1:ip.num_VRES)
            
                +
                
                sum( 
                    ip.conv.maintenance_costs[n,i,e]*
                    (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                    + 
                    ip.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
                for e in 1:ip.num_conv) 
                 
                for i in 1:ip.num_prod)
        
            - sum(

                ip.transm.maintenance_costs[n,m]*
                (ip.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m])
                +
                ip.transm.investment_costs[n,m]*0.5*l_plus[n,m]

            for m in 1:ip.num_nodes)
    
        for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )

    ## UPPER LOWEL CONSTRAINTS
    # Primal feasibility for the transmission expansion 1
    @constraint(Upper(bi_level_problem), [n in 1:ip.num_nodes, m in 1:ip.num_nodes],
        l_plus[n,m]/scaling_factor  - l_plus[m,n]/scaling_factor  == 0 
    )

    # Primal feasibility for the transmission expansion 2 (buget related)
        @constraint(Upper(bi_level_problem),
        sum(ip.transm.investment_costs[n,m] * 0.5 * l_plus[n,m] for n in 1:ip.num_nodes, m in 1:ip.num_nodes)/scaling_factor  - ip.transm.budget_limit/scaling_factor  <= 0 
    )
 

    ## LOWER LEVEL OBJECTIVE
    @objective(Lower(bi_level_problem), Max, 
            sum(
                sum( ip.scen_prob[s] * 

                    (ip.id_intercept[s,t,n]*q[s,t,n] 
                    - 0.5* ip.id_slope[s,t,n]/ip.time_periods[t]* q[s,t,n]^2


                    - (market == "perfect" ? 0 : 
                        (0.5* ip.id_slope[s,t,n]/ip.time_periods[t] * 
                        sum( 
                            (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2
                        for i in 1:ip.num_prod)
                        )
                    )
                    
                    - sum( 
                        (ip.conv.operational_costs[n,i,e] 
                        + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                    for i in 1:ip.num_prod, e in 1:ip.num_conv)

                    #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                    #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                    #- 1E-4 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                    )

                for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
                
                -sum( 

                    sum( 
                        ip.vres.maintenance_costs[n,i,r]*
                        (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                        + 
                        ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                    for r in 1:ip.num_VRES)
                    
                    +
                        
                    sum( 
                        ip.conv.maintenance_costs[n,i,e]*
                        (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                        + 
                        ip.conv.investment_costs[n,i,e]* g_conv_plus[n,i,e]
                    for e in 1:ip.num_conv) 
                         
                for i in 1:ip.num_prod)
            
            for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )

    
    ## LOWER LEVEL CONSTRAINTS 
    # Power balance
    @constraint(Lower(bi_level_problem), θ[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        q[s,t,n]/scaling_factor  - sum( 
                        sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
                         + 
                        sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
                        + g_hydro[s,t,n,i]
                        for i = 1:ip.num_prod)/scaling_factor 
                        + sum( f[s,t,n,m] for m in n+1:ip.num_nodes)/scaling_factor  - sum(f[s,t,m,n] for m in 1:n-1)/scaling_factor  #- sum( 0.98*f[s,t,m,n] for m in 1:n-1)
        == 0
    )

    # Transmission bounds 

    @constraint(Lower(bi_level_problem), β_f_1[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    @constraint(Lower(bi_level_problem), β_f_2[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        -f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    # Primal feasibility for the transmission 
    @constraint(Lower(bi_level_problem), λ_f[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        f[s,t,n,m]/scaling_factor  == 0 
    )


    # Generation budget limits
    @constraint(Lower(bi_level_problem), β_plus[i = 1:ip.num_prod],
        sum(ip.vres.investment_costs[n,i,r] * g_VRES_plus[n,i,r] for  r in 1:ip.num_VRES, n in 1:ip.num_nodes)/scaling_factor 
        + sum(ip.conv.investment_costs[n,i,e] * g_conv_plus[n,i,e] for  e in 1:ip.num_conv, n in 1:ip.num_nodes)/scaling_factor <= ip.gen_budget[i]/scaling_factor
    )


    # Conventional generation bounds
    @constraint(Lower(bi_level_problem), β_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t]*(ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    # VRES generation bounds 
    @constraint(Lower(bi_level_problem), β_VRES[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        g_VRES[s,t,n,i,r]/scaling_factor  - round(ip.time_periods[t]*ip.vres.availability_factor[s,t,n,r], digits = 4)*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r])/scaling_factor  <= 0
    )


    # Hydro power generation bounds 
    @constraint(Lower(bi_level_problem), [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod],
        g_hydro[s,t,n,i]/scaling_factor - (ip.hydro[s,t,n,i])/scaling_factor <= 0
    )

    # Maximum ramp-up rate for conventional units
    @constraint(Lower(bi_level_problem), β_up_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    # Maximum ramp-down rate for conventional units
    @constraint(Lower(bi_level_problem), β_down_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )
    
    
    return bi_level_problem
end

function lower_level_problem_generation_new(lines_investments_fixed::NTuple{10, Int64}, ip::initial_parameters, market::String)    
    # Defining single-level model
    lower_level_problem = Model(() -> Gurobi.Optimizer())
    #set_optimizer_attribute(lower_level_problem, "Threads", 1)
    set_optimizer_attribute(lower_level_problem, "OutputFlag", 0)

    ## UPPER LEVEL VARIABLES

    # Transmission capacity expansion realted variable
    l_plus = Array{Float64}(undef, ip.num_nodes, ip.num_nodes)
    for i in 1:length(lines_investments_fixed)
        l_plus[transmission_line_index_map[i]] = lines_investments_fixed[i]
        l_plus[transmission_line_index_map[i][2], transmission_line_index_map[i][1]] = lines_investments_fixed[i]
    end
    #print(l_plus)
    ## LOWER LEVEL VARIABLES 

    # Conventional generation related variable
    @variable(lower_level_problem, g_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # VRES generation related variable
    @variable(lower_level_problem, g_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)
    
    # hydro power generation related variable
    @variable(lower_level_problem, g_hydro[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod] >= 0)

    # Consumption realted variable
    @variable(lower_level_problem, q[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes]>= 0)

    # Energy transmission realted variable
    @variable(lower_level_problem, f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes] ) #>= 0)

    # Conventional energy capacity expansion realted variable
    @variable(lower_level_problem, g_conv_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv]>=0)

    # VRES capacity expansion realted variable
    @variable(lower_level_problem, g_VRES_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES]>=0)  

    ## LOWER LEVEL OBJECTIVE
    @objective(lower_level_problem, Max, 
            sum(
                sum( ip.scen_prob[s] * 

                    (ip.id_intercept[s,t,n]*q[s,t,n] 
                    - 0.5* ip.id_slope[s,t,n]/ip.time_periods[t]* q[s,t,n]^2


                    - (market == "perfect" ? 0 : 
                        (0.5* ip.id_slope[s,t,n]/ip.time_periods[t] * 
                        sum( 
                            (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2
                        for i in 1:ip.num_prod)
                        )
                    )
                    
                    - sum( 
                        (ip.conv.operational_costs[n,i,e] 
                        + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                    for i in 1:ip.num_prod, e in 1:ip.num_conv)

                    #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                    #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                    #- 1E-4 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                    )

                for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
                
                -sum( 

                    sum( 
                        ip.vres.maintenance_costs[n,i,r]*
                        (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                        + 
                        ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                    for r in 1:ip.num_VRES)
                    
                    +
                        
                    sum( 
                        ip.conv.maintenance_costs[n,i,e]*
                        (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                        + 
                        ip.conv.investment_costs[n,i,e]* g_conv_plus[n,i,e]
                    for e in 1:ip.num_conv) 
                         
                for i in 1:ip.num_prod)
            
            for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )

    
    ## LOWER LEVEL CONSTRAINTS 
    # Power balance
    @constraint(lower_level_problem, θ[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        q[s,t,n]/scaling_factor  - sum( 
                        sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
                         + 
                        sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
                        + g_hydro[s,t,n,i]
                        for i = 1:ip.num_prod)/scaling_factor 
                        + sum( f[s,t,n,m] for m in n+1:ip.num_nodes)/scaling_factor  - sum(f[s,t,m,n] for m in 1:n-1)/scaling_factor  #- sum( 0.98*f[s,t,m,n] for m in 1:n-1)
        == 0
    )

    # Transmission bounds 

    @constraint(lower_level_problem, β_f_1[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    @constraint(lower_level_problem, β_f_2[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        -f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    # Primal feasibility for the transmission 
    @constraint(lower_level_problem, λ_f[s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        f[s,t,n,m]/scaling_factor  == 0 
    )


    # Generation budget limits
    @constraint(lower_level_problem, β_plus[i = 1:ip.num_prod],
        sum(ip.vres.investment_costs[n,i,r] * g_VRES_plus[n,i,r] for  r in 1:ip.num_VRES, n in 1:ip.num_nodes)/scaling_factor 
        + sum(ip.conv.investment_costs[n,i,e] * g_conv_plus[n,i,e] for  e in 1:ip.num_conv, n in 1:ip.num_nodes)/scaling_factor <= ip.gen_budget[i]/scaling_factor
    )


    # Conventional generation bounds
    @constraint(lower_level_problem, β_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t]*(ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    # VRES generation bounds 
    @constraint(lower_level_problem, β_VRES[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        g_VRES[s,t,n,i,r]/scaling_factor  - round(ip.time_periods[t]*ip.vres.availability_factor[s,t,n,r], digits = 4)*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r])/scaling_factor  <= 0
    )

    # Hydro power generation bounds 
    @constraint(lower_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod],
        g_hydro[s,t,n,i]/scaling_factor - (ip.hydro[s,t,n,i])/scaling_factor <= 0
    )

    # Maximum ramp-up rate for conventional units
    @constraint(lower_level_problem, β_up_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    # Maximum ramp-down rate for conventional units
    @constraint(lower_level_problem, β_down_conv[ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )
    
    
    return lower_level_problem, l_plus
end