src_link  = "/Users/nikitabelyak/Dropbox (Aalto)/TSEP/Nordic_case_study"
cd(src_link)
using Pkg
Pkg.activate(".")

Pkg.instantiate()
using CSV, DataFrames, JuMP, Gurobi, PyCall, Plots,  Statistics, XLSX, BilevelJuMP, Cbc, Dates, Suppressor, Glob
#np = pyimport("numpy")

# set unique envinronment for Gurobi
const GRB_ENV = Gurobi.Env()

# scaling factor for the monetary parameters 
scaling_factor = 1.0

# scaling factor for the objective coefficients 

obg_scaling_factor = 1.0

# maximum number of decimal digits when roudning 
#( minimum value should correspond to the biggest numebr of digits used in scaling factors)

#max_digits = max(1, Int(max(log10(scaling_factor), log10(obg_scaling_factor))))
global global_max_digits = 3


data_src_link  =  src_link * "/2022 data"


EXP_INPUT_TAX = [[0.0 8.0 3.0 4.0 2.0], [0.0 70.0 25.0 36.0 16.0]]
EXP_GEN_BUDGET = [1153846.0, 57592307.0]
EXP_TRANSM_BUDGET = [1153846.0, 57592307.0]
EXP_INC = [5, 20]

mutable struct EXP_INPUT_SET
    gen_budget::String
    transm_budget::String
    tax::String
    incentive::String
end

global input_set = EXP_INPUT_SET("large", "large", "large", "large")

include(src_link*"/src/utils/data_preprocessing_functions.jl")
include(src_link*"/src/types/parameters.jl")
include(src_link*"/src/utils/parameters_initialisation.jl")
include(src_link*"/src/utils/models_generation_new.jl")
include(src_link*"/src/utils/data_postprocessing_functions.jl")

single_level_problem = single_level_problem_generation(input_parameters)
optimize!(single_level_problem)
print_output(data_src_link*"/optimisation_results" , input_parameters, objective_value(single_level_problem), value.(single_level_problem[:l_plus]), value.(single_level_problem[:g_VRES_plus]) , value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:g_VRES]), value.(single_level_problem[:g_conv]), value.(single_level_problem[:f]), value.(single_level_problem[:q]), incentives, input_parameters.gen_budget, "single_level", scaling_factor)
result_count(single_level_problem)

value.(single_level_problem[:l_plus])

bi_level_problem_perfect = bi_level_problem_generation_new(input_parameters, "perfect")
optimize!(bi_level_problem_perfect)
print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_perfect), value.(bi_level_problem_perfect[:l_plus]), value.(bi_level_problem_perfect[:g_VRES_plus]) , value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:g_VRES]), value.(bi_level_problem_perfect[:g_conv]), value.(bi_level_problem_perfect[:f]), value.(bi_level_problem_perfect[:q]), incentives, input_parameters.gen_budget, "bi_level_perfect", scaling_factor)
result_count(bi_level_problem_perfect)


# Maximum capacity investment (MW)
maxci = 9000
# Minimum capacity investment (MW)
minci = 0
# range step (MW)
step = 3000


# collect all the poossible combinations of the investemnts for 10 transmission lines 

all_inv_combinations = collect(Iterators.product(minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci));
all_inv_combinations = reshape(all_inv_combinations, 1, length(all_inv_combinations)); 


# the function checks the wether the termination status
# had any errors and returns true if it has and fals otherwise
function termination_status_error_check(status_code::Int64)
    # if INFEASIBLE or INFEASIBLE_OR_UNBOUNDED
    return  ((status_code == 2) || (status_code == 6))
end

global transmission_line_index_map = Dict(1 => CartesianIndex(1,2), 2 => CartesianIndex(1,3), 3 => CartesianIndex(1,4), 4 => CartesianIndex(1,5), 5 => CartesianIndex(2,3),
                                    6 => CartesianIndex(2,4), 7 => CartesianIndex(2,5), 8 => CartesianIndex(3,4), 9 => CartesianIndex(3,5), 10 => CartesianIndex(4,5))
global best_obj_value = zeros(length(all_inv_combinations))

Threads.nthreads() = Sys.CPU_THREADS



df1 = DataFrame(set_up=String[], obj_value_cen = Float64[], VRES_share_cen = Float64[], Total_generation_cen = Float64[], obj_value_per = Float64[], VRES_share_per = Float64[], Total_generation_per = Float64[] )

global input_set = EXP_INPUT_SET("small", "small", "small", "small")

include(src_link*"/src/utils/data_preprocessing_functions.jl")
include(src_link*"/src/types/parameters.jl")
include(src_link*"/src/utils/parameters_initialisation.jl")
include(src_link*"/src/utils/models_generation_new.jl")
include(src_link*"/src/utils/data_postprocessing_functions.jl")

all_inv_combinations = (0, 0, 3000, 0, 0, 3000, 0, 9000, 0, 3000)



for combination = 1:1

    lower_level_problem_perfect, l_plus_value = lower_level_problem_generation_new(all_inv_combinations, input_parameters, "perfect")
    @show l_plus_value
    bi_level_problem_perfect = bi_level_problem_generation_new(input_parameters, "perfect")
    single_level_problem = single_level_problem_generation(input_parameters)

    #@show input_parameters.transm.budget_limit - sum( input_parameters.transm.investment_costs[n,m]*0.5*l_plus_value[n,m] for n in  1:input_parameters.num_nodes, m in 1:input_parameters.num_nodes)
    TRANSM_EXPENSES = sum( input_parameters.transm.maintenance_costs[n,m]*(input_parameters.transm.installed_capacities[n,m] + 0.5 * l_plus_value[n,m]) + input_parameters.transm.investment_costs[n,m]*0.5*l_plus_value[n,m] for n in  1:input_parameters.num_nodes, m in 1:input_parameters.num_nodes)

    for i = 1:size(l_plus_value,1)
        for j = 1:size(l_plus_value,2)
            fix(bi_level_problem_perfect[:l_plus][i,j], l_plus_value[i,j], force=true)
            fix(single_level_problem[:l_plus][i,j], l_plus_value[i,j], force=true)
        end
    end
    optimize!(lower_level_problem_perfect)
    optimize!(single_level_problem)
    print("GEN BUDGET: $(input_set.gen_budget), TRANSM BUDGET: $(input_set.transm_budget), TAX: $(input_set.tax), INCENTIVE: $(input_set.incentive):\n\n")
    ov_central, vres_share_central, total_generation_central = print_output(data_src_link*"/optimisation_results" , input_parameters, objective_value(single_level_problem), value.(single_level_problem[:l_plus]), value.(single_level_problem[:g_VRES_plus]) , value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:g_VRES]), value.(single_level_problem[:g_conv]), value.(single_level_problem[:f]), value.(single_level_problem[:q]), incentives, input_parameters.gen_budget, "single_level", scaling_factor)
    #print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_perfect), value.(bi_level_problem_perfect[:l_plus]), value.(bi_level_problem_perfect[:g_VRES_plus]) , value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:g_VRES]), value.(bi_level_problem_perfect[:g_conv]), value.(bi_level_problem_perfect[:f]), value.(bi_level_problem_perfect[:q]), incentives, input_parameters.gen_budget, "bi_level_perfect", scaling_factor)
    print("CENTRALISED PLANNER:\n")
    print("welafre: $(ov_central)\n")
    print("VRES share: $(vres_share_central)\n")
    print("total generation: $(total_generation_central)\n\n")
    
    ov_perfect, vres_share_perfect, total_generation_perfect = print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(lower_level_problem_perfect), l_plus_value, value.(lower_level_problem_perfect[:g_VRES_plus]) , value.(lower_level_problem_perfect[:g_conv_plus]), value.(lower_level_problem_perfect[:g_VRES]), value.(lower_level_problem_perfect[:g_conv]), value.(lower_level_problem_perfect[:f]), value.(lower_level_problem_perfect[:q]), incentives, input_parameters.gen_budget, "bi_level_perfect", scaling_factor)
    print("PERFECT COMPETITION:\n")
    print("welafre: $(ov_perfect - TRANSM_EXPENSES)\n")
    print("VRES share: $(vres_share_perfect)\n")
    print("total generation: $(total_generation_perfect)\n\n")

    push!(df1, ("GEN BUDGET: $(input_set.gen_budget), TRANSM BUDGET: $(input_set.transm_budget), TAX: $(input_set.tax), INCENTIVE: $(input_set.incentive):\n\n", (ov_central), vres_share_central, total_generation_central, (ov_perfect - TRANSM_EXPENSES), vres_share_perfect, total_generation_perfect  ) )
end
CSV.write(folder*"/combined_results_1.csv", df1)
#df = DataFrame(function_value = collect(best_obj_value[1:100]))
#CSV.write(data_src_link * "/optimisation_results/lower_level_objective_value_$(input_set.gen_budget)_gen_b_$(input_set.transm_budget)_gen_b_$(input_set.tax)_tax_$(input_set.incentive)_inc.csv", df)


#detecting the indeces at which the total investment overshoot the allocated budget (for both cases with big and small budget)
global bad_indexes_small_budget = []
global bad_indexes_bigg_budget = []

for l = 1:length(all_inv_combinations)
    if sum(input_parameters.transm.investment_costs[transmission_line_index_map[i]]* all_inv_combinations[l][i] for i = 1:10) > 1153846.0
        append!(bad_indexes_small_budget,l)
    end

    if sum(input_parameters.transm.investment_costs[transmission_line_index_map[i]]* all_inv_combinations[l][i] for i = 1:10) > 57592307.0
        append!(bad_indexes_bigg_budget,l)
    end
end

# going through all the data with the welafre values for each of the input parameters set
folder = data_src_link * "/optimisation_results/triton_results"
files = glob("*.csv", folder) 
df = DataFrame(set_up=String[], obj_value=Float64[], investments=String[] )

for file in files
    obj_values_triton = Array(DataFrame.(CSV.File.(file)))
    
    # avoiding those indexes where the transmission investment exceeded the budget by making seure they will never be max. 
    if file[end-34:end-30] == "large" 
        obj_values_triton[bad_indexes_bigg_budget] .= -1000
    elseif file[end-34:end-30] == "small"
        @show file[end-34:end-30]
        obj_values_triton[bad_indexes_small_budget] .= -1000
    end
    vl, index = findmax(obj_values_triton)

    push!(df, (file[end-46:end-42]*", "*file[end-34:end-30]*", "*file[end-22:end-18]*", "*file[end-12:end-8], vl, string(all_inv_combinations[index[1]]) ) )

end

CSV.write(folder*"/combined_results.csv", df)


