cd(@__DIR__)

using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Printf

 
# Helper functions
 
function sample_name_from_inputfile(path)
    fname = splitpath(path)[end]
    fname = replace(fname, ".csv" => "")
    fname = replace(fname, "_ratios" => "")
    return strip(fname)
end

# Ion naming convention, replace for appropriate convention
function locate_norm_oxide_column(df, ox)
    for suffix in ("_oxide_norm_wt%", "_oxide_norm_wtpct")
        name = "$(ox)$(suffix)"
        name in names(df) && return name
    end
    return nothing
end

function mean_bulk_from_dataframe(df, oxide_cols, oxides)
    means = Dict{String,Float64}()
    for ox in oxides
        means[ox] = mean(df[!, oxide_cols[ox]])
    end
    total = sum(values(means))
    if total > 0
        for ox in oxides
            means[ox] = means[ox] / total * 100
        end
    end
    for k in keys(means)
        if isnan(means[k])
            means[k] = 0.0
        end
    end
    return means
end

function bulk_from_bulkgeochem(sample_name, bulkfile, oxides)
    bulkdf = CSV.read(bulkfile, DataFrame)
    sample_col = names(bulkdf)[1]
    rows = bulkdf[bulkdf[!, sample_col] .== sample_name, :]
    nrow(rows) == 0 && error("No bulk geochem found for sample $sample_name")
    row = rows[1, :]
    bulk = Dict{String,Float64}()
    for ox in oxides
        bulk[ox] = 0.0
    end
    for ox in ["SiO2","Al2O3","Fe2O3","MnO","MgO","CaO","Na2O","K2O","TiO2","P2O5"]
        if ox in names(bulkdf)
            bulk[ox] = Float64(row[ox]) * 100.0
        end
    end
    return bulk
end

# Oxide classification, replace with appropriate oxides
classification_oxides = [
    "SiO2","TiO2","Al2O3","Cr2O3","Fe2O3","MnO",
    "MgO","CaO","Na2O","K2O","P2O5","ZrO2"
]

#Euclidean distance function
function distance_to_mineral(measured, target)
    tgt = [get(target, ox, 0.0) for ox in classification_oxides]
    return norm(measured .- tgt)
end

# Normalise a mineral-formula Dict to sum to 100
function normalize_to_100!(formula)
    s = sum(values(formula))
    if s > 0
        for k in keys(formula)
            formula[k] = formula[k] / s * 100
        end
    end
    return formula
end

# Spot classification function

function classify_spot(measured, mineral_formulas, distfn)
    scored = [(m, distfn(measured, f)) for (m, f) in mineral_formulas]
    sort!(scored, by = x -> x[2])
    return scored[1][1], scored[1][2], scored[2][1], scored[2][2]
end

# Confidence rating from top-1 / top-2 ratio
function confidence_label(d1, d2)
    d1 == 0 && return "high"
    r = d2 / d1
    r >= 2.0 && return "high"
    r >= 1.3 && return "med"
    return "low"
end

# Finds files, replace with appropriate file paths
 
files = readdir(@__DIR__)
input_files = filter(f ->
    occursin("_ratios", f) &&
    endswith(f, ".csv")     &&
    !occursin("STANDARDS", f),
    files
)
println("Found $(length(input_files)) ratio files.")

 
# Mineral identification loop
 
for input_file in input_files
    println("\n==============================")
    println("Processing $input_file")

    sample_name       = sample_name_from_inputfile(input_file)
    output_file       = "$(sample_name)_minerals.csv"                      # Replace with appropriate naming scheme, biotite spots can be removed, but replaced with any other different types of spots measured
    biotite_output    = "$(sample_name)_biotite_spots.csv"
    bulk_geochem_file = "bulk_geochem.csv"

    println("Reading $input_file …")
    df_all = CSV.read(input_file, DataFrame)

 
    # Split off biotite spots into a separate file
 
    df_biotite = DataFrame()
    df         = df_all

    if "Spot_ID" in names(df_all)
        is_biotite = occursin.(r"_biotite_"i, string.(df_all.Spot_ID))
        is_grid    = occursin.(r"_grid_\d+$"i, string.(df_all.Spot_ID))

        df_biotite = df_all[is_biotite, :]
        df         = df_all[is_grid, :]

        println("Rows: total=$(nrow(df_all))  grid=$(nrow(df))  biotite=$(nrow(df_biotite))")
    else
        println("⚠️ No Spot_ID column found; using all rows for classification.")
    end

 
    # Find  and clean oxides column 
 
    oxide_cols = Dict{String,String}()
    for ox in classification_oxides
        col = locate_norm_oxide_column(df, ox)
        if col === nothing
            cname = "$(ox)_oxide_norm_wt%"
            df[!, cname] = zeros(Float64, nrow(df))
            oxide_cols[ox] = cname
        else
            oxide_cols[ox] = col
        end
    end
 
    for col in values(oxide_cols)
        df[!, col] = coalesce.(df[!, col], 0.0)
        df[!, col] = ifelse.(isnan.(Float64.(df[!, col])), 0.0, Float64.(df[!, col]))
    end

 
    # Re-normalize selected oxides to 100%
 
    selected_cols = [oxide_cols[o] for o in classification_oxides]
    selected_sum  = sum(Matrix(select(df, selected_cols)), dims=2)[:,1]
    for col in selected_cols
        newcol = Symbol(replace(col, r"_oxide_norm_wt(%|pct)" => "_selected_norm_wt%"))
        df[!, newcol] = ifelse.(selected_sum .> 0,
                                df[!, col] ./ selected_sum .* 100,
                                0.0)
    end
    selected_norm_cols = [replace(c, r"_oxide_norm_wt(%|pct)" => "_selected_norm_wt%") for c in selected_cols]

 
    # Mineral formulas list - add minerals depending on desired identification, and change depending on rock type
 
    mineral_formulas = Dict{String,Dict{String,Float64}}(

    "biotite" => Dict(
    "SiO2"=>37.0,"TiO2"=>3.5,"Al2O3"=>14.0,"Cr2O3"=>0.0,
    "Fe2O3"=>16.0,"MnO"=>0.3,"MgO"=>15.0,"CaO"=>0.0,
    "Na2O"=>0.3,"K2O"=>9.5,"P2O5"=>0.0,"ZrO2"=>0.0),

    "muscovite" => Dict(
    "SiO2"=>45.0,"TiO2"=>0.0,"Al2O3"=>38.0,"Cr2O3"=>0.0,
    "Fe2O3"=>1.0,"MnO"=>0.0,"MgO"=>0.5,"CaO"=>0.0,
    "Na2O"=>0.5,"K2O"=>11.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "sericite" => Dict(
    "SiO2"=>48.0,"TiO2"=>0.2,"Al2O3"=>32.0,"Cr2O3"=>0.0,
    "Fe2O3"=>3.0,"MnO"=>0.0,"MgO"=>2.0,"CaO"=>0.2,
    "Na2O"=>0.3,"K2O"=>10.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "calcite" => Dict(
    "SiO2"=>0.0,"TiO2"=>0.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.0,"MnO"=>0.0,"MgO"=>0.5,"CaO"=>56.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "plagioclase-ca" => Dict(
    "SiO2"=>52.0,"TiO2"=>0.0,"Al2O3"=>30.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.3,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>12.0,
    "Na2O"=>4.5,"K2O"=>0.3,"P2O5"=>0.0,"ZrO2"=>0.0),

    "plagioclase-na" => Dict(
    "SiO2"=>62.0,"TiO2"=>0.0,"Al2O3"=>24.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.2,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>4.0,
    "Na2O"=>9.0,"K2O"=>0.5,"P2O5"=>0.0,"ZrO2"=>0.0),

    "potassium feldspar" => Dict(
    "SiO2"=>65.0,"TiO2"=>0.0,"Al2O3"=>18.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.2,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>0.2,
    "Na2O"=>1.0,"K2O"=>15.5,"P2O5"=>0.0,"ZrO2"=>0.0),

    "chlorite" => Dict(
    "SiO2"=>29.0,"TiO2"=>0.1,"Al2O3"=>20.0,"Cr2O3"=>0.0,
    "Fe2O3"=>22.0,"MnO"=>0.5,"MgO"=>22.0,"CaO"=>0.2,
    "Na2O"=>0.0,"K2O"=>0.1,"P2O5"=>0.0,"ZrO2"=>0.0),

    "hornblende" => Dict(
    "SiO2"=>46.0,"TiO2"=>1.4,"Al2O3"=>9.0,"Cr2O3"=>0.0,
    "Fe2O3"=>13.0,"MnO"=>0.4,"MgO"=>14.0,"CaO"=>11.5,
    "Na2O"=>1.5,"K2O"=>0.8,"P2O5"=>0.0,"ZrO2"=>0.0),

    "actinolite" => Dict(
    "SiO2"=>55.0,"TiO2"=>0.1,"Al2O3"=>2.0,"Cr2O3"=>0.0,
    "Fe2O3"=>7.0,"MnO"=>0.3,"MgO"=>21.0,"CaO"=>13.0,
    "Na2O"=>0.2,"K2O"=>0.1,"P2O5"=>0.0,"ZrO2"=>0.0),

    "diopside" => Dict(
    "SiO2"=>53.0,"TiO2"=>0.2,"Al2O3"=>1.0,"Cr2O3"=>0.0,
    "Fe2O3"=>7.0,"MnO"=>0.2,"MgO"=>17.0,"CaO"=>24.0,
    "Na2O"=>0.2,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "quartz" => Dict(
    "SiO2"=>99.5,"TiO2"=>0.0,"Al2O3"=>0.2,"Cr2O3"=>0.0,
    "Fe2O3"=>0.0,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>0.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.3),

    "sphene" => Dict(
    "SiO2"=>30.5,"TiO2"=>37.0,"Al2O3"=>1.5,"Cr2O3"=>0.0,
    "Fe2O3"=>1.5,"MnO"=>0.1,"MgO"=>0.0,"CaO"=>28.5,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>1.0),

    "rutile" => Dict(
    "SiO2"=>0.0,"TiO2"=>99.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.5,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>0.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "apatite" => Dict(
    "SiO2"=>0.5,"TiO2"=>0.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.2,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>54.0,
    "Na2O"=>0.3,"K2O"=>0.0,"P2O5"=>42.0,"ZrO2"=>0.0),

    "zircon" => Dict(
    "SiO2"=>33.0,"TiO2"=>0.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.0,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>0.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>66.0),

    "iron oxide" => Dict(
    "SiO2"=>0.2,"TiO2"=>1.5,"Al2O3"=>0.3,"Cr2O3"=>0.0,
    "Fe2O3"=>97.0,"MnO"=>0.3,"MgO"=>0.2,"CaO"=>0.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "ilmenite" => Dict(
    "SiO2"=>0.0,"TiO2"=>52.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>46.0,"MnO"=>1.5,"MgO"=>0.5,"CaO"=>0.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "orthopyroxene" => Dict(
    "SiO2"=>54.0,"TiO2"=>0.2,"Al2O3"=>2.0,"Cr2O3"=>0.0,
    "Fe2O3"=>20.0,"MnO"=>0.5,"MgO"=>22.0,"CaO"=>1.5,
    "Na2O"=>0.1,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    "anhydrite" => Dict(
    "SiO2"=>0.0,"TiO2"=>0.0,"Al2O3"=>0.0,"Cr2O3"=>0.0,
    "Fe2O3"=>0.0,"MnO"=>0.0,"MgO"=>0.0,"CaO"=>41.0,
    "Na2O"=>0.0,"K2O"=>0.0,"P2O5"=>0.0,"ZrO2"=>0.0),

    )
    # Renormalise mineral formulas to 100 wt%
    for (name, f) in mineral_formulas
        normalize_to_100!(f)
    end

 
    # Glass/matrix identification - uses bulk formula here, but change depending on desired matrix identification scheme
 
    println("Detected sample: '$sample_name'")                         #For samples without known bulk, finds mean values
    if startswith(sample_name, "EM")
        println("✅ Using EM mean bulk")                                                                 
        bulk_mean = mean_bulk_from_dataframe(df, oxide_cols, classification_oxides)
        mineral_formulas["glass/bulk"] = bulk_mean
    else
        println("✅ Using bulk_geochem.csv")
        mineral_formulas["glass/bulk"] = bulk_from_bulkgeochem(sample_name, bulk_geochem_file, classification_oxides)
    end
    
    excluded_from_euclidean = ["calcite", "anhydrite"]
    mineral_formulas_euclidean = Dict(k => v for (k, v) in mineral_formulas
                                      if !(k in excluded_from_euclidean))

 
    # Euclidean identification, without calcite or anhydrite
 
    n = nrow(df)
    df.Top1_Mineral              = Vector{String}(undef, n)
    df.Top1_Distance             = Vector{Float64}(undef, n)
    df.Top2_Mineral              = Vector{String}(undef, n)
    df.Top2_Distance             = Vector{Float64}(undef, n)
    df.Confidence                = Vector{String}(undef, n)

    for i in 1:n
        measured = Float64[df[i, col] for col in selected_norm_cols]
        b1, d1, b2, d2 = classify_spot(measured, mineral_formulas_euclidean, distance_to_mineral)
        df.Top1_Mineral[i]   = b1
        df.Top1_Distance[i]  = d1
        df.Top2_Mineral[i]   = b2
        df.Top2_Distance[i]  = d2
        df.Confidence[i]     = confidence_label(d1, d2)
    end

    # Save results for grid points
 
    println("Writing $output_file …")
    CSV.write(output_file, df)

 
    # Save results for other (biotite) points
 
    if nrow(df_biotite) > 0

        biotite_oxide_cols = Dict{String,String}()
        for ox in classification_oxides
            col = locate_norm_oxide_column(df_biotite, ox)
            if col === nothing
                cname = "$(ox)_oxide_norm_wt%"
                df_biotite[!, cname] = zeros(Float64, nrow(df_biotite))
                biotite_oxide_cols[ox] = cname
            else
                biotite_oxide_cols[ox] = col
            end
        end

        for col in values(biotite_oxide_cols)
            df_biotite[!, col] = coalesce.(df_biotite[!, col], 0.0)
            df_biotite[!, col] = ifelse.(isnan.(Float64.(df_biotite[!, col])), 0.0, Float64.(df_biotite[!, col]))
        end

        biotite_selected_cols = [biotite_oxide_cols[o] for o in classification_oxides]
        biotite_selected_sum  = sum(Matrix(select(df_biotite, biotite_selected_cols)), dims=2)[:,1]
        for col in biotite_selected_cols
            newcol = Symbol(replace(col, r"_oxide_norm_wt(%|pct)" => "_selected_norm_wt%"))
            df_biotite[!, newcol] = ifelse.(biotite_selected_sum .> 0,
                                            df_biotite[!, col] ./ biotite_selected_sum .* 100,
                                            0.0)
        end
        biotite_norm_cols = [replace(c, r"_oxide_norm_wt(%|pct)" => "_selected_norm_wt%") for c in biotite_selected_cols]

        n_b = nrow(df_biotite)
        df_biotite.Top1_Mineral                     = Vector{String}(undef, n_b)
        df_biotite.Top1_Distance                    = Vector{Float64}(undef, n_b)
        df_biotite.Top2_Mineral                     = Vector{String}(undef, n_b)
        df_biotite.Top2_Distance                    = Vector{Float64}(undef, n_b)
        df_biotite.Confidence                       = Vector{String}(undef, n_b)

        #Euclidean identification, without calcite or anhydrite
        for i in 1:n_b
            measured = Float64[df_biotite[i, col] for col in biotite_norm_cols]

            b1, d1, b2, d2 = classify_spot(measured, mineral_formulas_euclidean, distance_to_mineral)
            df_biotite.Top1_Mineral[i]  = b1
            df_biotite.Top1_Distance[i] = d1
            df_biotite.Top2_Mineral[i]  = b2
            df_biotite.Top2_Distance[i] = d2
            df_biotite.Confidence[i]    = confidence_label(d1, d2)

        end

    end

end  

println("\n✅ ALL FILES CLASSIFIED")