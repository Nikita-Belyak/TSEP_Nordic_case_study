import Pkg
Pkg.activate(".")
Pkg.add("TimeSeriesClustering")
using TimeSeriesClustering
data_path=normpath(joinpath(dirname(@__FILE__),"tutorial"))
your_data_2 = load_timeseries_data(data_path; T=24, years=[2018])