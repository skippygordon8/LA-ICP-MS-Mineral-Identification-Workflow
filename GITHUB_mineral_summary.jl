using CSV
using DataFrames
using Printf
using Glob

cd(@__DIR__)

#Change for apppropriate file names
input_folder          = @__DIR__
output_file_euclidean = "mineral_summary_euclidean.csv"
output_file_biotite   = "mineral_summary_biotite.csv"

 
# Maps identified minerals to common/useful names, change for appropriate names
 
function mineral_to_category(mineral::AbstractString)
    if mineral == "glass/bulk"
        return "matrix"
    elseif mineral == "mode"
        return "mode"
    elseif mineral == "sericite" || mineral == "muscovite"
        return "muscovite"
    elseif mineral == "chlorite"
        return "chlorite"
    elseif mineral == "plagioclase-ca" || mineral == "plagioclase-na"
        return "plagioclase"
    elseif mineral == "potassium feldspar"
        return "orthoclase"
    elseif mineral == "quartz"
        return "quartz"
    elseif mineral == "actinolite" || mineral == "hornblende"
        return "amphibole"
    elseif mineral == "biotite"
        return "biotite"
    elseif mineral == "sphene"
        return "sphene"
    elseif mineral == "apatite"
        return "apatite"
    elseif mineral == "orthopyroxene"
        return "orthopyroxene"
    elseif mineral == "iron oxide" || mineral == "ilmenite" || mineral == "rutile"
        return "opaques"
    elseif mineral == "diopside"
        return "clinopyroxene"
    elseif mineral == "epidote" || mineral == "epidote-pistacite"
        return "epidote"
    elseif mineral == "calcite"
        return "calcite"
    elseif mineral == "zircon"
        return "zircon"
    elseif mineral == "anhydrite"
        return "anhydrite"
    else
        return "matrix"
    end
end

 
# Find sample name from file path, change as appropriate
 
function sample_name_from_filename(path::AbstractString)
    fname = splitpath(path)[end]
    fname = replace(fname, "_biotite_spots.csv" => "")
    fname = replace(fname, "_minerals.csv" => "")
    return fname
end

 
# Finds files in directory, change as appropriate
 
all_files     = glob("*.csv", input_folder)
grid_files    = filter(f -> occursin("_minerals", lowercase(f)) && !occursin("biotite", lowercase(f)) && !occursin("summary", lowercase(f)), all_files)
biotite_files = filter(f -> occursin("_biotite_spots", lowercase(f)), all_files)

isempty(grid_files) && error("No grid *_minerals.csv files found in $input_folder")

 
# Prepare output columns, change as appropriate
 
col_schema = (
    sample        = String[],
    unit          = String[],
    matrix        = Int[],
    mode          = Int[],
    plagioclase   = Int[],
    orthoclase    = Int[],
    quartz        = Int[],
    amphibole     = Int[],
    biotite       = Int[],
    chlorite      = Int[],
    muscovite     = Int[],
    sphene        = Int[],
    apatite       = Int[],
    opaques       = Int[],
    clinopyroxene = Int[],
    orthopyroxene = Int[],
    epidote       = Int[],
    zircon        = Int[],
    calcite       = Int[],
    anhydrite     = Int[],
    total         = Int[],
)

summary_euclidean = DataFrame(col_schema)
summary_biotite   = DataFrame(col_schema)

 
# Function for adding mineral counts
 
function count_and_push!(summary, df, mineral_col, sample)
    counts = Dict(
        "matrix"        => 0,
        "mode"          => 0,
        "plagioclase"   => 0,
        "orthoclase"    => 0,
        "quartz"        => 0,
        "amphibole"     => 0,
        "biotite"       => 0,
        "chlorite"      => 0,
        "muscovite"     => 0,
        "sphene"        => 0,
        "apatite"       => 0,
        "opaques"       => 0,
        "clinopyroxene" => 0,
        "orthopyroxene" => 0,
        "epidote"       => 0,
        "zircon"        => 0,
        "calcite"       => 0,
        "anhydrite"     => 0,
    )
    for mineral in df[!, mineral_col]
        cat = mineral_to_category(string(mineral))
        counts[cat] += 1
    end
    total = sum(values(counts))
    push!(summary, (
        sample, "",
        counts["matrix"],
        counts["mode"],
        counts["plagioclase"],
        counts["orthoclase"],
        counts["quartz"],
        counts["amphibole"],
        counts["biotite"],
        counts["chlorite"],
        counts["muscovite"],
        counts["sphene"],
        counts["apatite"],
        counts["opaques"],
        counts["clinopyroxene"],
        counts["orthopyroxene"],
        counts["epidote"],
        counts["zircon"],
        counts["calcite"],
        counts["anhydrite"],
        total
    ))
end

 
# Process point files
 
for file in grid_files
    println("Reading $file ...")
    df = CSV.read(file, DataFrame)

    if "Spot_ID" in names(df)
        df = filter(row -> occursin("grid", lowercase(string(row.Spot_ID))), df)
    else
        @printf("⚠️  No Spot_ID column in %s — skipping\n", file)
        continue
    end

    if nrow(df) == 0
        @printf("⚠️  No grid points in %s — skipping\n", file)
        continue
    end

    if !("Top1_Mineral" in names(df)) || !("Mineral_Diagnostic_Weighted" in names(df))
        @printf("⚠️  Skipping %s because mineral columns are missing\n", file)
        continue
    end

    sample = sample_name_from_filename(file)
    count_and_push!(summary_euclidean, df, "Top1_Mineral",                sample)
    count_and_push!(summary_weighted,  df, "Mineral_Diagnostic_Weighted", sample)
end

 
# Process biotite point files
for file in biotite_files
    println("Reading $file ...")
    df = CSV.read(file, DataFrame)

    if nrow(df) == 0
        @printf("⚠️  No rows in %s — skipping\n", file)
        continue
    end

    if !("Mineral_Diagnostic_Weighted" in names(df))
        @printf("⚠️  Skipping %s because Mineral_Diagnostic_Weighted column is missing\n", file)
        continue
    end

    sample = sample_name_from_filename(file)
    count_and_push!(summary_biotite, df, "Mineral_Diagnostic_Weighted", sample)
end

 
# Write output file
 
println("Writing $output_file_euclidean ...")
CSV.write(output_file_euclidean, summary_euclidean)

println("Writing $output_file_weighted ...")
CSV.write(output_file_weighted, summary_weighted)

if nrow(summary_biotite) > 0
    println("Writing $output_file_biotite ...")
    CSV.write(output_file_biotite, summary_biotite)
else
    println("(No biotite files found — skipping biotite summary)")
end

println("✅ Done")