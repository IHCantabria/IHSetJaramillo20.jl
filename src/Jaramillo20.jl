"""
#######################################
# Jaramillo et al. (2020) model       #
# E: wave energy                      #
# dt: time step                       #
# a: model's parameter                #
# b: model's parameter                #
# kero: erosion coefficient           #
# kacr: accretion coefficient         #
# Yi: initial position                #
# vlt: longshore velocity             #
#######################################
"""


function run_Jaramillo20()

    println("Loading libraries...")
    wrkDir = pwd()

    println("Loading datasets...")

    wavF = dats*"wav.nc"
    parF = dats*"par.nc"
    slF = dats*"sl.nc"
    configF = dats*"config.nc"

    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt = configF["dt"][:][1]
    
    yi = configF["yi"][:][1]

    a = parF["a"][:][1]
    b = parF["b"][:][1]
    cacr = parF["cacr"][:][1]
    cero = parF["cero"][:][1]
    

    
    brk, angBati, depth, Hberm, D50 = configF["brk"][:][1], configF["angBati"][:][1], configF["depth"][:][1], configF["Hberm"][:][1], configF["D50"][:][1]

    if brk == 1
        
        Hs, Tp, θ_w = collect(skipmissing(wavF["Hs"][:])), collect(skipmissing(wavF["Tp"][:])), collect(skipmissing(wavF["Dir"][:]))

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, Tp, Hs, depthb = wavF["Hb"][:], wavF["Tp"][:], wavF["Hs"][:], wavF["hb"][:]
    end
    
    close(wavF)
    close(ensF)
    close(configF)
    close(parF)

    println("Datasets closed...")

    Hs = convert(Array{Float64},Hs)
    Tp = convert(Array{Float64},Tp)

    a = convert(Array{Float64},a)
    b = convert(Array{Float64},b)
    cacr = convert(Array{Float64},kacr)
    cero = convert(Array{Float64},kero)

    Yi = convert(Array{Float64},yi)

    ########## START HERE #############

    E = Hb .^ 2

    println("Starting Jaramillo et al. (2020) - Shoreline Rotation Model...")


    Y_t = Jaramillo20(E, dt, a, b, cacr, cero, Yi)

    println("\n\n****************Finished****************\n\n")

    return Y_t
    
end

function Jaramillo20(E, dt, a, b, cacr, cero, Yini, vlt)
    
    Seq = (E .- b) ./ a

    Y = zeros(length(E))
    
    Y[1] = Yini

    for i in eachindex(E[2:end])
        if Y[i] < Seq[i+1]
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cacr *(E[i+1] ^ 0.5)*dt))+Seq[i+1] + vlt*dt
        else
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cero *(E[i+1] ^ 0.5)*dt))+Seq[i+1] + vlt*dt
        end
    end
    
    return Y
end


function cal_Jaramillo20()

    """
    cal_Jaramillo20
    
    Function to calibrate and run the Jaramillo et al. (2020) Shoreline Evolution Model.
    
    This function reads input datasets, performs calibration, and writes the results to an output NetCDF file.
    
    Note: The function internally uses the Yates09 function for shoreline evolution.
    
    """

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    parF = dats*"par.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]

    calYi = ncread(configF, "calYi")[1]
    
    calVlt = ncread(configF, "calVlt")[1]

    brk, angBati, depth, D50 = ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1], ncread(configF, "D50")[1]

    MetObj = ncread(configF, "MetObj")[1]

    if brk == 1
        
        Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

        auxAng, auxDepth = zeros(length(Hs)), zeros(length(Hs))
        auxAng .= angBati
        auxDepth .= depth
 
        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    elseif brk == 2
        Hb, Tp, Hs = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Hs")
        θ_b = ncread(wavF, "Dir")
        depthb = zeros(length(Hb)) .+ depth
    else
        Hb, Tp, Hs, depthb = ncread(wavF, "Hb"), ncread(wavF, "Tp"), ncread(wavF, "Hs"), ncread(wavF, "hb")
    end

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, hs, depthb = t_wav[ii], Hb[ii], Tp[ii], Hs[ii], depthb[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)

    ########## START HERE #############

    E = Hb .^ 2

    println("Starting Jaramillo et al. (2020) - Shoreline Evolution Model...")

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()

    if calYi == 0 && calVlt == 0
        vlt = ncread(configF, "vlt")[1]
        function Calibra_4r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Y_obs[1], vlt)
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_4r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), Y_obs[1], vlt)

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)


        println("\n\n****************Finished****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")


        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([vlt], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")

    elseif calYi == 1 && calVlt == 0 
        vlt = ncread(configF, "vlt")[1]
        function Calibra_5r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Χ[5], vlt)
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (0.25*minimum(Y_obs), 2*maximum(Y_obs))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_5r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), popr[5], vlt)

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([vlt], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMS",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE_flagP="*string(i))
        nccreate(output, "MSS_flagP="*string(i),
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS_flagP="*string(i))

    elseif calYi == 1 && calVlt == 1
        function Calibra_6r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Χ[5], Χ[6])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (0.25*minimum(Y_obs), 2*maximum(Y_obs)),
                    (-100/(365.25*24), 100/(365.25*24))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_6r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_6r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_6r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), popr[5], popr[6])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([popr[6]], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")

    elseif calYi == 0 && calVlt == 1
        function Calibra_7r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]),Y_obs[1], Χ[5])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (-100/(365.25*24), 100/(365.25*24))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_7r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_7r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_7r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), Y_obs[1], popr[5])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([popr[5]], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")

    end


    println("\n\n****************Finished****************\n\n")

end


function calVal_Jaramillo20()

    """
    calVal_Jaramillo20
    
    Function to calibrate, in a given calibration period, and validate the Jaramillo et al. (2020) Shoreline Evolution Model.
    
    This function reads input datasets, performs calibration, and writes the results to an output NetCDF file.
    
    Note: The function internally uses the Yates09 function for shoreline evolution.
    
    """

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    parF = dats*"par.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]

    calYi = ncread(configF, "calYi")[1]
    
    calVlt = ncread(configF, "calVlt")[1]

    brk, angBati, depth, D50 = ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1], ncread(configF, "D50")[1]

    MetObj = ncread(configF, "MetObj")[1]

    if brk == 1
        
        Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

        auxAng, auxDepth = zeros(length(Hs)), zeros(length(Hs))
        auxAng .= angBati
        auxDepth .= depth
 
        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    elseif brk == 2
        Hb, Tp, Hs = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Hs")
        θ_b = ncread(wavF, "Dir")
        depthb = zeros(length(Hb)) .+ depth
    else
        Hb, Tp, Hs, depthb = ncread(wavF, "Hb"), ncread(wavF, "Tp"), ncread(wavF, "Hs"), ncread(wavF, "hb")
    end

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")

    YYs, MMs, DDs, HHs = ncread(configF, "Ys"), ncread(configF, "Ms"), ncread(configF, "Ds"), ncread(configF, "hs")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)
    t_split = DateTime.(YYs, MMs, DDs, HHs)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, hs, depthb = t_wav[ii], Hb[ii], Tp[ii], Hs[ii], depthb[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)

    idx_obs_split = Int(argmin(abs.(t_obs .- t_split)))

    idx_split = Int(argmin(abs.(t_wav .- t_split)))

    t_obs_validation, Y_obs_validation, idx_obs_validation = t_obs[idx_obs_split+1:end], Y_obs[idx_obs_split+1:end], idx_obs[idx_obs_split+1:end]

    t_obs, Y_obs, idx_obs = t_obs[1:idx_obs_split], Y_obs[1:idx_obs_split], idx_obs[1:idx_obs_split]

    ########## START HERE #############

    E = Hb .^ 2

    println("Starting Jaramillo et al. (2020) - Shoreline Evolution Model...")

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()

    if calYi == 0 && calVlt == 0
        vlt = ncread(configF, "vlt")[1]
        function Calibra_4r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Y_obs[1], vlt)
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_4r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), Y_obs[1], vlt)

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)
        
        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_v",
                    "dim", length(YYs),
                    atts = yearv_atts)
        ncwrite(YYs, output, "year_v")
        nccreate(output, "month_v",
                    "dim", length(MMs),
                    atts = monthv_atts)
        ncwrite(MMs, output, "month_v")
        nccreate(output, "day_v",
                    "dim", length(DDs),
                    atts = dayv_atts)
        ncwrite(DDs, output, "day_v")
        nccreate(output, "hour_v",
                    "dim", length(HHs),
                    atts = hourv_atts)
        ncwrite(HHs, output, "hour_v")
        
        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")

        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")


        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([vlt], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")


    elseif calYi == 1 && calVlt == 0 
        vlt = ncread(configF, "vlt")[1]
        function Calibra_5r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Χ[5], vlt)
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (0.25*minimum(Y_obs), 2*maximum(Y_obs))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_5r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), popr[5], vlt)

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_v",
                    "dim", length(YYs),
                    atts = yearv_atts)
        ncwrite(YYs, output, "year_v")
        nccreate(output, "month_v",
                    "dim", length(MMs),
                    atts = monthv_atts)
        ncwrite(MMs, output, "month_v")
        nccreate(output, "day_v",
                    "dim", length(DDs),
                    atts = dayv_atts)
        ncwrite(DDs, output, "day_v")
        nccreate(output, "hour_v",
                    "dim", length(HHs),
                    atts = hourv_atts)
        ncwrite(HHs, output, "hour_v")

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([vlt], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMS",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS=")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")

    elseif calYi == 1 && calVlt == 1
        function Calibra_6r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Χ[5], Χ[6])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (0.25*minimum(Y_obs), 2*maximum(Y_obs)),
                    (-100/(365.25*24), 100/(365.25*24))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_6r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_6r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_6r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 6,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 100000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), popr[5], popr[6])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_v",
                    "dim", length(YYs),
                    atts = yearv_atts)
        ncwrite(YYs, output, "year_v")
        nccreate(output, "month_v",
                    "dim", length(MMs),
                    atts = monthv_atts)
        ncwrite(MMs, output, "month_v")
        nccreate(output, "day_v",
                    "dim", length(DDs),
                    atts = dayv_atts)
        ncwrite(DDs, output, "day_v")
        nccreate(output, "hour_v",
                    "dim", length(HHs),
                    atts = hourv_atts)
        ncwrite(HHs, output, "hour_v")

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([popr[6]], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")

    elseif calYi == 0 && calVlt == 1
        function Calibra_7r(Χ)
            Ymd = Jaramillo20(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]),Y_obs[1], Χ[5])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (-100/(365.25*24), 100/(365.25*24))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_7r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_7r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_7r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 50000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 1000000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo20(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), Y_obs[1], popr[5])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA20.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_v",
                    "dim", length(YYs),
                    atts = yearv_atts)
        ncwrite(YYs, output, "year_v")
        nccreate(output, "month_v",
                    "dim", length(MMs),
                    atts = monthv_atts)
        ncwrite(MMs, output, "month_v")
        nccreate(output, "day_v",
                    "dim", length(DDs),
                    atts = dayv_atts)
        ncwrite(DDs, output, "day_v")
        nccreate(output, "hour_v",
                    "dim", length(HHs),
                    atts = hourv_atts)
        ncwrite(HHs, output, "hour_v")

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        vlt_atts = Dict("units" => "m/s",
            "long_name" => "Longshore velocity",
            "standard_name" => "vlt")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "vlt",
                    "len", 1,
                    atts = vlt_atts)
        ncwrite([popr[5]], output, "vlt")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")

    end


    println("\n\n****************Finished****************\n\n")

end