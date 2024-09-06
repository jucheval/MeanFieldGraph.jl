using DiscreteHawkes
using DataFrames
using Distributions
using ProgressLogging
using StatsPlots
using CSV
using TableMetadataTools
using Logging

begin   # Auxiliary functions
    function abserror(data::DiscreteHawkesData, r₊, Δ, truevalues)
        if isnan(Δ)
            Δ = max(floor(Int,log(length(data))), 1)
        else
            Δ = Int(Δ)
        end
        m̂, v̂, ŵ = estimators(data, Δ)
    
        μλphat = collect(DiscreteHawkes.Φ(m̂, v̂, ŵ, r₊))
        mvwhat = [m̂, v̂, ŵ]
        return abs.([mvwhat; μλphat] .- truevalues)
    end
    
    function abserror(modelconnec::ConnectivityMatrix, Z, r₊, truevalues)
        m̂, v̂, ŵ = mvw_inf(modelconnec, Z)
    
        μλphat = collect(DiscreteHawkes.Φ(m̂, v̂, ŵ, r₊))
        mvwhat = [m̂, v̂, ŵ]
        return abs.([mvwhat; μλphat] .- truevalues)
    end
    
    function complete!(df, E, Paramvec, tvec, plow, pup)
        n1, n2, n3, n4 = size(E)
        for idparam in 1:n4
            for idt in 1:n2
                q1 = quantile(E[1,idt,:,idparam], [plow,.5,pup])
                q2 = quantile(E[2,idt,:,idparam], [plow,.5,pup])
                q3 = quantile(E[3,idt,:,idparam], [plow,.5,pup])
                q4 = quantile(E[4,idt,:,idparam], [plow,.5,pup])
                q5 = quantile(E[5,idt,:,idparam], [plow,.5,pup])
                q6 = quantile(E[6,idt,:,idparam], [plow,.5,pup])
                Paramvec[idparam] == NaN ? param = 0 : param = Paramvec[idparam]
                push!(df, (idparam, param, tvec[idt], q1[1], q1[2], q1[3],
                                                                q2[1], q2[2], q2[3],
                                                                q3[1], q3[2], q3[3],
                                                                q4[1], q4[2], q4[3],
                                                                q5[1], q5[2], q5[3],
                                                                q6[1], q6[2], q6[3]))
            end
        end
    end
    
    function complete_inf!(df_inf, E_inf, Paramvec, plow, pup)
        n = size(E_inf,3)
        for idparam in 1:n
            q1 = quantile(E_inf[1,:,idparam], [plow,.5,pup])
            q2 = quantile(E_inf[2,:,idparam], [plow,.5,pup])
            q3 = quantile(E_inf[3,:,idparam], [plow,.5,pup])
            q4 = quantile(E_inf[4,:,idparam], [plow,.5,pup])
            q5 = quantile(E_inf[5,:,idparam], [plow,.5,pup])
            q6 = quantile(E_inf[6,:,idparam], [plow,.5,pup])
            Paramvec[idparam] == NaN ? param = 0 : param = Paramvec[idparam]
            push!(df, (idparam, param, q1[1], q1[2], q1[3],
                                                        q2[1], q2[2], q2[3],
                                                        q3[1], q3[2], q3[3],
                                                        q4[1], q4[2], q4[3],
                                                        q5[1], q5[2], q5[3],
                                                        q6[1], q6[2], q6[3]))
        end
    end
    
    function metadatacomplete!(df::DataFrame, Paramsymbol::Symbol, type, N, r₊, β, λ, p, Nsimu, plow, pup, Δ)
        metadata!(df, "Caption", "Performance of the DiscreteHawkes estimators")
        Paramsymbol == :N || metadata!(df, "N", N)
        Paramsymbol == :r₊ || metadata!(df, "r₊", r₊)
        Paramsymbol == :β || metadata!(df, "β", β)
        Paramsymbol == :λ || metadata!(df, "λ", λ)
        Paramsymbol == :p || metadata!(df, "p", p)
        metadata!(df, "Number of simulations", Nsimu)
        metadata!(df, "Quantiles proba.", [plow, .5, pup])
        metadata!(df, "Type of error", type)
        if "T" in names(df)
            if Paramsymbol != :Δ 
                if isnan(Δ)
                    metadata!(df, "Δ", "Chosen as the logarithm of T")
                else
                    metadata!(df, "Δ", Δ)
                end
            end
            label!(df, :T, "Time horizon used for the estimation")
        end
        metadata!(df, "Varying parameter", String(Paramsymbol))
        label!(df, :color, "Color used for plotting")
        label!(df, :parameter, "Value of the parameter "*string(Paramsymbol))
        for targetstring in ["m","v","w","μ","λ","p"]
            label!(df, Symbol(targetstring*"low"), string(100*plow)*"% quantile of the estimation error of "*targetstring)
            label!(df, Symbol(targetstring*"50"), string(100*0.5)*"% quantile of the estimation error of "*targetstring)
            label!(df, Symbol(targetstring*"up"), string(100*pup)*"% quantile of the estimation error of "*targetstring)
        end
    end
    
    function estimatorstable(model::DiscreteHawkesModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int})::Tuple{DataFrame, DataFrame}
        lenΔvec = length(Δvec)
        df = DataFrame(idsimu = Int[], T = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Vector{Float64}[])
        df_inf = DataFrame(idsimu = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[])
        
        Z = DiscreteHawkes.N2Z(N, r₊)
        
        @progress "simulation" for idsimu in 1:Nsimu
            θ = rand(DiscreteHawkes.ErdosRenyiGraph(N, model.p))
            modelconnec = DiscreteHawkes.ConnectivityMatrix(model,θ)
            push!(df_inf, (idsimu, mvw_inf(modelconnec, Z)...))

            data = rand(modelconnec, Z, maximum(tvec))
            for t in tvec
                tmpdata = data[1:t]
                push!(df, (idsimu, t, estimators(tmpdata, Δvec)...))
            end
        end
        return df, df_inf
    end

    function errortables(Paramsymbol::Symbol, type::String, Paramvec, Nsimu, tvec, Δ)
        type in ["absolute","relative"] || error("Argument type must be equal to absolute or relative")
        E = Array{Float64}(undef, 6, length(tvec), Nsimu, length(Paramvec))
        E_inf = Array{Float64}(undef, 6, Nsimu, length(Paramvec))
        for idparam in eachindex(Paramvec)
            if Paramsymbol == :N
                N = Paramvec[idparam]; global r₊, β, λ, p
            elseif Paramsymbol == :r₊
                r₊ = Paramvec[idparam]; global N, β, λ, p
            elseif Paramsymbol == :β
                β = Paramvec[idparam]; global r₊, N, λ, p
            elseif Paramsymbol == :λ
                λ = Paramvec[idparam]; global r₊, β, N, p
            elseif Paramsymbol == :p
                p = Paramvec[idparam]; global r₊, β, λ, N
            end
            model = DiscreteHawkesModel(β*λ,λ,p)
            m, v, w = mvw(model, r₊)
            truevalues = [m, v, w, model.μ, model.λ, model.p]
            Z = DiscreteHawkes.N2Z(N, r₊)
            if type == "relative"
                boundstruevalues = [[0,0,0,0,0,0] [1,1/4,10,1,1,1]]
                scaling = minimum(abs.(truevalues .- boundstruevalues),dims=2)
            else
                scaling = 1
            end
            @progress "sim and estim, i="*string(idparam)*" on "*string(length(Paramvec)) for idsimu in 1:Nsimu
                θ = rand(DiscreteHawkes.ErdosRenyiGraph(N, model.p))
                modelconnec = DiscreteHawkes.ConnectivityMatrix(model,θ)
                E_inf[:, idsimu, idparam] = abserror(modelconnec, Z, r₊, truevalues)
    
                data = rand(modelconnec, Z, maximum(tvec))
                for idt in eachindex(tvec)
                    if Paramsymbol == :Δ
                        E[:, idt, idsimu, idparam] = abserror(data[1:tvec[idt]], r₊, Paramvec[idparam], truevalues) ./ scaling 
                    else
                        E[:, idt, idsimu, idparam] = abserror(data[1:tvec[idt]], r₊, Δ, truevalues) ./ scaling
                    end
                end
            end
        end
    
        Paramsymbol in (:N, :Δ) ? Typeparam = Int : Typeparam = Float64
        df = DataFrame(color = Int[], parameter = Typeparam[], T = Int[], mlow = Float64[], m50 = Float64[], mup = Float64[],
                                                            vlow = Float64[], v50 = Float64[], vup = Float64[],
                                                            wlow = Float64[], w50 = Float64[], wup = Float64[],
                                                            μlow = Float64[], μ50 = Float64[], μup = Float64[],
                                                            λlow = Float64[], λ50 = Float64[], λup = Float64[],
                                                            plow = Float64[], p50 = Float64[], pup = Float64[])
        df_inf = DataFrame(color = Int[], parameter = Typeparam[], mlow = Float64[], m50 = Float64[], mup = Float64[],
                                                vlow = Float64[], v50 = Float64[], vup = Float64[],
                                                wlow = Float64[], w50 = Float64[], wup = Float64[],
                                                μlow = Float64[], μ50 = Float64[], μup = Float64[],
                                                λlow = Float64[], λ50 = Float64[], λup = Float64[],
                                                plow = Float64[], p50 = Float64[], pup = Float64[])
        complete!(df, E, Paramvec, tvec, plow, pup)
        metadatacomplete!(df, Paramsymbol, type, N, r₊, β, λ, p, Nsimu, plow, pup, Δ)
        complete_inf!(df_inf, E_inf, Paramvec, plow, pup)
        metadatacomplete!(df_inf, Paramsymbol, type, N, r₊, β, λ, p, Nsimu, plow, pup, Δ)
    
        return df, df_inf
    end
end;


# Default values
N = 500
r₊ = .5
β = .5
λ = .5
p = .5
T = Int(1e3)
Nsimu = Int(1e3)

Δ = 1
type = "absolute"

length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))
plow = .25
pup = .75

disable_logging(LogLevel(-1001))  # Enables debug info

begin   # Simulation 
    # Change the values and the symbol of the parameter
    Paramvec = [1, 2, 5, NaN]
    df, df_inf = errortables(:Δ, type, Paramvec, Nsimu, tvec, Δ)

    paramstring = metadata(df, "Varying parameter")
    # begin   # Save
    #     CSV.write("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*".csv", df)
    #     open("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*".toml", "w") do io
    #         print(io, meta2toml(df))
    #     end
    
    #     CSV.write("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.csv", df_inf)
    #     open("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.toml", "w") do io
    #         print(io, meta2toml(df_inf))
    #     end
    # end;
end;

begin   # Load
    # Change the param string and Delta value
    paramstring = "λ"
    Δ = 1
    
    df = CSV.read("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*".csv", DataFrame)
    open("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*".toml") do io
        toml2meta!(df, io)
    end
    df_inf = CSV.read("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.csv", DataFrame)
    open("../discrete-hawkes/code/data/error_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.toml") do io
        toml2meta!(df_inf, io)
    end
end;

begin   # 1 plot for each estimator
    for targetstring in ["m","v","w","μ","λ","p"]
        paramstring = metadata(df, "Varying parameter")
        plot(df.T, df[:,targetstring*"50"], group=df.parameter, color=df.color)
        #plot!(df.T, df[:,targetstring*"low"], group=df.parameter, color=df.color, label=false, alpha = .5, style = :dash)
        #plot!(df.T, df[:,targetstring*"up"], group=df.parameter, color=df.color, label=false, alpha = .5, style = :dash)
        scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring*"50"], group=df_inf.parameter, color=df_inf.color, marker = :hline, label=false, markerstrokewidth = 5)
        #scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring*"low"], group=df_inf.parameter, color=df_inf.color, marker = :dtriangle, label=false)
        #scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring*"up"], group=df_inf.parameter, color=df_inf.color, marker = :utriangle, label=false)
        xlabel!("T")
        ylims!(0, min(maximum(df[:,targetstring*"up"]),.3))
        ylabel!(metadata(df)["Type of error"]*" error")
        title!("Estimation error for "*targetstring*" as "*paramstring*" varies, Δ="*string(Δ))
        savefig("../discrete-hawkes/figures/outputs/error_"*targetstring*"_vary_"*paramstring*"_delta_"*string(Δ)*".pdf")
    end
end


begin # Simulation with default values
    df, df_inf = errortables(:N, type, [N], Nsimu, tvec, Δ)

    begin   # Save
        CSV.write("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*".csv", df)
        open("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*".toml", "w") do io
            print(io, meta2toml(df))
        end
    
        CSV.write("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*"_inf.csv", df_inf)
        open("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*"_inf.toml", "w") do io
            print(io, meta2toml(df_inf))
        end
    end;
end

begin   # 1 plot for all estimators
    # Load
    df = CSV.read("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*".csv", DataFrame)
    open("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*".toml") do io
        toml2meta!(df, io)
    end
    df_inf = CSV.read("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*"_inf.csv", DataFrame)
    open("../discrete-hawkes/code/data/error_default_values_delta_"*string(Δ)*"_inf.toml") do io
        toml2meta!(df_inf, io)
    end
    plot(yaxis=:log, legend_columns=3)
    for (id,targetstring) in enumerate(["m","v","w","μ","λ","p"])
        plot!(df.T, df[:,targetstring*"50"], color=id, label=targetstring)
        scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring*"50"], color=id, marker = :o, label=false)
    end
    xlabel!("T")
    ylims!(1e-4, 1)
    ylabel!(metadata(df)["Type of error"]*" error")
    title!("Estimation error for all parameters, Δ="*string(Δ))
    savefig("../discrete-hawkes/figures/outputs/error_default_values_delta_"*string(Δ)*".pdf")
end