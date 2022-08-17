src_link  = "/scratch/work/belyakn1/TSEP_nordic_case_study"
cd(src_link)
using Pkg
Pkg.activate(".")

#Pkg.update()

#ENV["GUROBI_HOME"] = "/share/apps/anaconda-ci/fgci-centos7-generic/software/anaconda/2020-01-tf2/5a34a04a"
#ENV["GUROBI_HOME"] = "/share/apps/anaconda-ci/fgci-centos7-anaconda/software/anaconda/2022-02/73aef705"
#ENV["GRB_LICENSE_FILE"]="/scratch/work/belyakn1/TSEP_nordic_case_study/gurobi.lic"

Pkg.add("Gurobi")
Pkg.build("Gurobi")

Pkg.instantiate()
using CSV, DataFrames, JuMP, Gurobi, PyCall,  Statistics, XLSX, BilevelJuMP, Dates, Suppressor, Base.Threads
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
EXP_TRANSM_BUDGET = [1153846.0, 57692307.0]
EXP_INC = [5, 20]

mutable struct EXP_INPUT_SET
    tax::String
    gen_budget::String
    transm_budget::String
    incentive::String
end

global input_set = EXP_INPUT_SET("small", "large", "small", "small")


include(src_link*"/src/utils/data_preprocessing_functions.jl")
include(src_link*"/src/types/parameters.jl")
include(src_link*"/src/utils/parameters_initialisation.jl")
include(src_link*"/src/utils/models_generation_new.jl")
include(src_link*"/src/utils/data_postprocessing_functions.jl")


#single_level_problem = single_level_problem_generation(input_parameters)
#optimize!(single_level_problem)
#print_output(data_src_link*"/optimisation_results" , input_parameters, objective_value(single_level_problem), value.(single_level_problem[:l_plus]), value.(single_level_problem[:g_VRES_plus]) , value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:g_VRES]), value.(single_level_problem[:g_conv]), value.(single_level_problem[:f]), value.(single_level_problem[:q]), incentives, input_parameters.gen_budget, "single_level", scaling_factor)
#result_count(single_level_problem)

#bi_level_problem_perfect = bi_level_problem_generation_new(input_parameters, "perfect")
#optimize!(bi_level_problem_perfect)
#print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_perfect), value.(bi_level_problem_perfect[:l_plus]), value.(bi_level_problem_perfect[:g_VRES_plus]) , value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:g_VRES]), value.(bi_level_problem_perfect[:g_conv]), value.(bi_level_problem_perfect[:f]), value.(bi_level_problem_perfect[:q]), incentives, input_parameters.gen_budget, "bi_level_perfect", scaling_factor)
#result_count(bi_level_problem_perfect)

# collect all the poossible combinations of the investemnts for 10 transmission lines 

# Maximum capacity investment (MW)
maxci = 9000
# Minimum capacity investment (MW)
minci = 0
# range step (MW)
step = 3000

all_inv_combinations = collect(Iterators.product(minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci, minci:step:maxci));
all_inv_combinations = reshape(all_inv_combinations, 1, length(all_inv_combinations)); 

global transmission_line_index_map = Dict(1 => CartesianIndex(1,2), 2 => CartesianIndex(1,3), 3 => CartesianIndex(1,4), 4 => CartesianIndex(1,5), 5 => CartesianIndex(2,3),
                                    6 => CartesianIndex(2,4), 7 => CartesianIndex(2,5), 8 => CartesianIndex(3,4), 9 => CartesianIndex(3,5), 10 => CartesianIndex(4,5))
global best_obj_value = zeros(length(all_inv_combinations))

# the function checks the wether the termination status
# had any errors and returns true if it has and fals otherwise
function termination_status_error_check(status_code::Int64)
    # if INFEASIBLE or INFEASIBLE_OR_UNBOUNDED
    return  ((status_code == 2) || (status_code == 6))
end


time_start = time()
@sync Base.Threads.@threads for combination = 1:length(all_inv_combinations)
    lower_level_problem_perfect, l_plus = lower_level_problem_generation_new(all_inv_combinations[combination], input_parameters, "perfect")
    if sum(l_plus)/2 <= input_parameters.transm.budget_limit 
        TRANSM_EXPENSES = sum( input_parameters.transm.maintenance_costs[n,m]*(input_parameters.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m]) + input_parameters.transm.investment_costs[n,m]*0.5*l_plus[n,m] for n in  1:input_parameters.num_nodes, m in 1:input_parameters.num_nodes)
        status = optimize!(lower_level_problem_perfect)
        if !termination_status_error_check(Int(termination_status(lower_level_problem_perfect)))
            best_obj_value[combination] = objective_value(lower_level_problem_perfect) - TRANSM_EXPENSES
        end
    end
    #@show time() - time_start
    #@show combination
end
@show time() - time_start
df = DataFrame(function_value = collect(best_obj_value[1:end]))
CSV.write(data_src_link * "/optimisation_results/lower_level_objective_value_$(input_set.gen_budget)_gen_b_$(input_set.transm_budget)_gen_b_$(input_set.tax)_tax_$(input_set.incentive)_inc.csv", df)


