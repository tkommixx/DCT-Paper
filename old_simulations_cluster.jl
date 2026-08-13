using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using Dates
using DelimitedFiles

#? cluster parameters and packages
using ClusterManagers
ENV["JULIA_WORKER_TIMEOUT"] = "180"
addprocs(ClusterManagers.SlurmManager(250), N=8, topology=:master_worker, exeflags="--project=Project.toml"; W="300")

## load the packages by old_covid19abm

#using old_covid19abm

#addprocs(10, exeflags="--project=.")


#@everywhere using old_covid19abm

# Taiye (2025.06.01): Return addprocs when connecting to the cluster
# Taiye (Use half cluster for testing): addprocs(ClusterManagers.SlurmManager(500), N=16, topology=:master_worker, exeflags = "--project=.")
# Taiye (2025.06.21): addprocs(ClusterManagers.SlurmManager(250), N=8, topology=:master_worker, exeflags = "--project=.")

@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("old_covid19abm.jl")
@everywhere const cv=old_covid19abm


function run(myp::cv.old_ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    #dump(myp)
   
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors     

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)

    # Taiye (2025.06.12): We are not considering workplaces.
    # working = vcat([cdr[i].work for i = 1:nsims]...)
   
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
   # ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    #ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
 #   ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
  #  ag7 = vcat([cdr[i].g7 for i = 1:nsims]...)

    # Taiye (2025.08.06):
    use_1 = vcat([cdr[i].u1 for i = 1:nsims]...)
    use_2 = vcat([cdr[i].u2 for i = 1:nsims]...)

    # Taiye (2025.06.12): We are not considering workplaces.
    # mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6,"ag7" => ag7, "working"=>working)
    mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "use_1" => use_1, "use_2" => use_2)

    #println(mydfs["all"]) # Taiye (2025.06.21): Testing
    #mydfs = Dict("all" => allag, "working"=>working, "kids"=>kids)
    #mydfs = Dict("all" => allag)
    
    #c1 = Symbol.((:LAT, :PRE, :MILD, :INF, :HOS, :ICU, :DED,:LAT2, :PRE2, :MILD2, :INF2, :HOS2, :ICU2, :DED2,:LAT3, :PRE3, :MILD3, :INF3, :HOS3, :ICU3, :DED3), :_INC)
    
    # Taiye:
    # c1 = Symbol.((:LAT, :MILD, :INF, :HOS, :ICU, :DED), :_INC)
    c1 = Symbol.((:LAT, :INF, :DED), :_INC)

    #c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2,:LAT3, :HOS3, :ICU3, :DED3), :_PREV)
    
    #c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        #for c in vcat(c1..., c2...)
        for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end

    
    writedlm(string(folderprefix,"/R01.dat"),[cdr[i].R0 for i=1:nsims])
    writedlm(string(folderprefix,"/year_of_death.dat"),hcat([cdr[i].vector_dead for i=1:nsims]...))
    # Taiye (2025.06.06): writedlm(string(folderprefix,"/npcr.dat"),hcat([cdr[i].npcr for i=1:nsims]...))
    writedlm(string(folderprefix,"/nra.dat"),hcat([cdr[i].nra for i=1:nsims]...))
    # Taiye (2025.06.06): writedlm(string(folderprefix,"/nleft.dat"),hcat([cdr[i].nleft for i=1:nsims]...))
    writedlm(string(folderprefix,"/totalisog.dat"),vcat([cdr[i].giso for i=1:nsims]))
    #writedlm(string(folderprefix,"/totalisow.dat"),vcat([cdr[i].wiso for i=1:nsims]))

    # Taiye (2025.08.06):
    #? Thomas: Here, you do not need the vcat, because you are returning one single number. Just like the R0.
    writedlm(string(folderprefix,"/totalquar.dat"),[cdr[i].quar_tot for i=1:nsims])

    # Taiye (2025.12.12):
    writedlm(string(folderprefix,"/min_age.dat"),[cdr[i].min_age for i=1:nsims])
    writedlm(string(folderprefix,"/app_age.dat"),[cdr[i].app_age for i=1:nsims])
    writedlm(string(folderprefix,"/ret_age.dat"),[cdr[i].ret_age for i=1:nsims])

   # writedlm(string(folderprefix,"/lat_inc.dat"),hcat([cdr[i].lat_inc' for i=1:nsims]...))
#    writedlm(string(folderprefix,"/asymp_inc.dat"),hcat([cdr[i].asymp_inc' for i=1:nsims]...))
 #   writedlm(string(folderprefix,"/pre_inc.dat"),hcat([cdr[i].pre_inc' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/inf_inc.dat"),hcat([cdr[i].inf_inc' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/inf_app_inc.dat"),hcat([cdr[i].inf_app_inc' for i=1:nsims]...))
    #writedlm(string(folderprefix,"/asymp_prev.dat"),hcat([cdr[i].asymp_prev' for i=1:nsims]...))
#    writedlm(string(folderprefix,"/pre_prev.dat"),hcat([cdr[i].pre_prev' for i=1:nsims]...))
 #   writedlm(string(folderprefix,"/inf_prev.dat"),hcat([cdr[i].inf_prev' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/inf_app_prev.dat"),hcat([cdr[i].inf_app_prev' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/rec_prev.dat"),hcat([cdr[i].rec_prev' for i=1:nsims]...))
    #writedlm(string(folderprefix,"/ded_prev.dat"),hcat([cdr[i].ded_prev' for i=1:nsims]...))
    
    # Taiye (2025.10.23):
 #   writedlm(string(folderprefix,"/entry_lat.dat"),hcat([cdr[i].entry_lat' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/entry_asymp.dat"),hcat([cdr[i].entry_asymp' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/entry_pre.dat"),hcat([cdr[i].entry_pre' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/entry_inf.dat"),hcat([cdr[i].entry_inf' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/entry_inf_app.dat"),hcat([cdr[i].entry_inf_app' for i=1:nsims]...))

  #  writedlm(string(folderprefix,"/dur_asymp.dat"),hcat([cdr[i].dur_asymp' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/dur_pre.dat"),hcat([cdr[i].dur_pre' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/dur_inf.dat"),hcat([cdr[i].dur_inf' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/dur_inf_app.dat"),hcat([cdr[i].dur_inf_app' for i=1:nsims]...))
  #  writedlm(string(folderprefix,"/dur_rec.dat"),hcat([cdr[i].dur_rec' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/dur_ded.dat"),hcat([cdr[i].dur_ded' for i=1:nsims]...))

  #  writedlm(string(folderprefix,"/lat_inc.dat"),hcat([cdr[i].lat_inc' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/asymp_inc.dat"),hcat([cdr[i].asymp_inc' for i=1:nsims]...))
    #writedlm(string(folderprefix,"/pre_inc.dat"),hcat([cdr[i].pre_inc' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/inf_inc.dat"),hcat([cdr[i].inf_inc' for i=1:nsims]...))
    #writedlm(string(folderprefix,"/inf_app_inc.dat"),hcat([cdr[i].inf_app_inc' for i=1:nsims]...))

    #writedlm(string(folderprefix,"/lat_prev.dat"),hcat([cdr[i].lat_prev' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/asymp_prev.dat"),hcat([cdr[i].asymp_prev' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/pre_prev.dat"),hcat([cdr[i].pre_prev' for i=1:nsims]...))
   # writedlm(string(folderprefix,"/inf_prev.dat"),hcat([cdr[i].inf_prev' for i=1:nsims]...))
    #writedlm(string(folderprefix,"/inf_app_prev.dat"),hcat([cdr[i].inf_app_prev' for i=1:nsims]...))


    return mydfs
end


function create_folder(ip::cv.old_ModelParameters,province="ontario")
    
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    #Taiye (2025.05.30): main_folder = "/data/thomas-covid/testing_canada"
    main_folder = "/data/Taiye"
    #main_folder = "."
    
    # Taiye (2025.05.27): secondaryfolder = string(main_folder,"/fmild_$(ip.fmild)_fwork_$(ip.fwork)") ##  
    #secondaryfolder = string(main_folder,"/fmild_$(ip.fmild)") # fwork is not found in old_covid19abm.jl
    
    # Taiye (2025.05.27): RF = string(secondaryfolder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_herd_immu_","$(ip.herd)","_idx_$(ip.file_index)_$(province)_strain_$(ip.strain)_scen_$(ip.scenariotest)_test_$(ip.test_ra)_eb_$(ip.extra_booster)_size_$(ip.size_threshold)") ## 
    # Taiye (2025.06.23): RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_idx_$(ip.file_index)_$(province)") 
    
    # Taiye (2025.07.01): Adding date to folder.
    # RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_idx_$(ip.file_index)_$(province)","_cov_$(ip.app_coverage)") 
    RF = string(main_folder,"/no_mandatory_symp_iso_notif_$(ip.not_swit)","_06_14_n_tests_$(ip.n_tests)","_results_prob_","$(replace(string(round(ip.β,digits=5)), "." => "_"))","_cov_$(round(ip.app_coverage,digits=2))","_comp_$(ip.comp_bool)","_int_$(ip.time_until_testing)","_iso_$(ip.iso_con)") 

    if !Base.Filesystem.isdir(main_folder)
        Base.Filesystem.mkpath(main_folder)
    end

    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end


# change coverage of app_coverage
# time testing
# Taiye (2025.05.27):
# function run_param_scen_cal(b::Float64,province::String="ontario",h_i::Int64 = 0,ic1::Int64=1,strains::Int64 = 1,index::Int64 = 0,scen::Int64 = 0,tra::Int64 = 0,eb::Int64 = 0,wpt::Int64 = 100,mt::Int64=300,test_time::Int64 = 1,test_dur::Int64=112,mildcomp::Float64 = 1.0,workcomp::Float64 = 1.0,dayst::Vector{Int64} = [1;4],trans_omicron::Float64 = 1.0,immu_omicron::Float64 = 0.0,rc=[1.0],dc=[1],nsims::Int64=500)
function run_param_scen_cal(b::Float64,province::String="ontario",ic1::Int64=1,index::Int64 = 0,test_time::Int64=0,test_dur::Int64=0,mt::Int64=300,nsims::Int64=500,ps::Int64=100000, app_cov::Float64=0.3, test_cap::Int64=1, i_con::Int64=0,notif::Bool=true,t_sens::Int64=1, ob::Bool=true, ims::Int64=1)
      
    @everywhere ip = cv.old_ModelParameters(β=$b,
    initialinf = $ic1,
    file_index = $index,
    modeltime=$mt, prov = Symbol($province),
    start_testing = $test_time,
    test_for = $test_dur,
    popsize = $ps,
    app_coverage = $app_cov,
    num_sims = $nsims,
    n_tests = $test_cap,
    iso_con = $i_con,
    not_swit = $notif,
    test_sens = $t_sens,
    comp_bool = $ob,
    time_until_testing = $ims
    )

    folder = create_folder(ip,province)

    run(ip,nsims,folder)
    #run(ip,4,folder)
   
end

bta = 0.12
run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,5,1,true,0,true,1)
run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,6,1,true,0,true,1) 
run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,7,1,true,0,true,1) 
run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,10,1,true,0,true,1) 

#run_param_scen_cal(bta,"ontario",1,0,1,365,365,10,10000,1.0,1,1,true,0,true,1) 


#run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,1,1,true,0,true,1)
#run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,1,1,true,0,true,2)
#run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,2,1,true,0,true,1)
#run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,1.0,2,1,true,0,true,2)


#for i = 1:10

    #run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.1*i,1,1,true,0,true,1)
   # run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.1*i,1,1,true,0,true,2)
  #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.1*i,2,1,true,0,true,1)
 #   run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.1*i,2,1,true,0,true,2)
    
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,false,1) # 1 test, 50% compliance, same day test
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,false,2) # 1 test, 50% compliance, next day test
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,false,3) # 1 test, 50% compliance, int = 3
#   run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,false,4) # 1 test, 50% compliance, int = 4      
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,false,1) # 2 tests, 50% compliance, same day test
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,false,2) # 2 tests, 50% compliance, next day test
   #run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,false,3) # 2 tests, 50% compliance, int = 3
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,false,4) # 2 tests, 50% compliance, int = 4
              
   
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.2*i,1,1,true,0,true,1) # 1 test, 100% compliance, same day test
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.2*i,1,1,true,0,true,2) # 1 test, 100% compliance, next day test
#   run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,true,3) # 1 test, 100% compliance, 3 days later
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,1,1,true,0,true,4) # 1 test, 100% compliance, 4 days later
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.2*i,2,1,true,0,true,1) # 2 tests, 100% compliance, same day test
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,1000,10000,0.2*i,2,1,true,0,true,2) # 2 tests, 100% compliance, next day test
 #  run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,true,3) # 2 tests, 100% compliance, 3 days later
  # run_param_scen_cal(bta,"ontario",1,0,1,365,365,500,10000,0.2*i,2,1,true,0,true,4) # 2 tests, 100% compliance, 4 days later
#   # end
#end