using CSV
using DataFrames
using Statistics
using Glob
cd(@__DIR__)

input_folder = @__DIR__
output_suffix = "_ratios.csv" #Change to desired suffix

 
# Search folder for files, change for appropriate naming scheme

files = glob("*_processed.csv", input_folder)
println("Found $(length(files)) files")

 
# Element list, change for appropriate elements
 
elements = [
    "Li", "Na", "Mg", "Al", "Si", "P", "K", "Ca", "Sc", "Ti",
    "V", "Cr", "Mn", "Fe", "Ni", "Cu", "Zn", "Rb", "Sr", "Y",
    "Zr", "Nb", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta"
]

 
# Selected oxides, change for appropriate major oxides  (ox_name, ox_mw, element_aw, n_cations)
 
oxide_map = Dict(
    "Si" => ("SiO2",     60.0843,  28.085,  1),
    "Ti" => ("TiO2",     79.866,   47.867,  1),
    "Al" => ("Al2O3",   101.9613,  26.982,  2),
    "Cr" => ("Cr2O3",   151.9904,  51.996,  2),
    "Fe" => ("Fe2O3",   159.687,   55.845,  2),
    "Mn" => ("MnO",      70.937,   54.938,  1),
    "Mg" => ("MgO",      40.304,   24.305,  1),
    "Ca" => ("CaO",      56.077,   40.078,  1),
    "Na" => ("Na2O",     61.9789,  22.99,   2),
    "K"  => ("K2O",      94.196,   39.098,  2),
    "P"  => ("P2O5",    141.9445,  30.974,  2),
    "Zr" => ("ZrO2",    123.222,   91.224,  1),       #Not major, but for the case a zircon got zapped
)

 
# Sample name extraction to prepare each sample for individual csvs, change to match appropriate sample naming scheme

function extract_sample(spot_id)
    s = String(spot_id)

    m = match(r"^(G_NIST\d+)", s)
    m !== nothing && return m.captures[1]

    m = match(r"^(HM22-\d+)", s)
    m !== nothing && return m.captures[1]

    m = match(r"^(EM24-[A-Za-z]+)", s)
    m !== nothing && return m.captures[1]

    return "UNKNOWN"
end

 
# Counts to ratios loop
 
for file in files
    println("\nProcessing $file ...")

    df = CSV.read(file, DataFrame)

 
    # Find element columns (strict regex: <Elem><mass>_ppm)
 
    element_cols = Dict{String,Symbol}()
    for elem in elements
        pattern = Regex("^" * elem * raw"\d+_ppm$")
        for col in names(df)
            if occursin(pattern, String(col))
                element_cols[elem] = Symbol(col)
                break
            end
        end
    end

    if isempty(element_cols)
        println("⚠️ No matching element columns — skipping")
        continue
    end

    println("Matched ", length(element_cols), " elements: ",
            join(sort(collect(keys(element_cols))), ", "))

 
    # CLEAN in temporary columns (preserve originals)
 
    clean_cols = Dict{String,Symbol}()
    for (elem, col) in element_cols
        newcol = Symbol(elem * "_ppm_clean")
        df[!, newcol] = max.(coalesce.(df[!, col], 0.0), 0.0)
        clean_cols[elem] = newcol
    end

 
    # Ion wt%  (true wt% = ppm / 10 000)
 
    for (elem, col) in clean_cols
        outcol = Symbol(elem * "_ion_wt%")
        df[!, outcol] = df[!, col] ./ 10_000.0
    end

 
    # Oxide conversion
 
    oxide_ppm_cols = String[]

    for (elem, (ox, ox_mw, el_aw, ncat)) in oxide_map
        if haskey(clean_cols, elem)
            col  = clean_cols[elem]
            conv = ox_mw / (el_aw * ncat)

            ppm_col = Symbol(ox * "_oxide_ppm")
            df[!, ppm_col] = df[!, col] .* conv

            wt_col = Symbol(ox * "_oxide_wt%")
            df[!, wt_col] = df[!, ppm_col] ./ 10_000.0

            push!(oxide_ppm_cols, String(ppm_col))
        end
    end

 
    # Normalize Oxides to 100 wt%
 
    if !isempty(oxide_ppm_cols)
        total_oxide = sum(Matrix(select(df, Symbol.(oxide_ppm_cols))), dims=2)[:, 1]

        for col in oxide_ppm_cols
            ox = replace(col, "_oxide_ppm" => "")
            outcol = Symbol(ox * "_oxide_norm_wtpct")

            df[!, outcol] = ifelse.(total_oxide .> 0,
                                    df[!, Symbol(col)] ./ total_oxide .* 100,
                                    0.0)
        end
    end

 
    # Drop temporary _clean columns
 
    select!(df, Not([v for v in values(clean_cols)]))

 
    # Writes csv for each sample 
 
    if !("Spot_ID" in names(df))
        @warn "No Spot_ID column found in $file — skipping (cannot split by sample)"
        continue
    end

    df.Sample = extract_sample.(df.Spot_ID)

    grouped = groupby(df, :Sample)
    println("Found $(length(grouped)) sample group(s) in $file:")

    for sub in grouped
        sample = sub.Sample[1]
        nrows  = nrow(sub)

        # Skip NIST standards and unknowns
        if startswith(sample, "G_NIST")
            println("  • $sample  ($nrows rows)  → SKIPPED (standard)")
            continue
        end
        if sample == "UNKNOWN"
            println("  • UNKNOWN     ($nrows rows)  → SKIPPED (unparseable Spot_ID)")
            continue
        end

        outfile = "$(sample)_ratios.csv"

        if isfile(outfile)
            println("  • $sample  ($nrows rows)  → appending to $outfile")
            existing = CSV.read(outfile, DataFrame)
            combined = vcat(existing, sub; cols = :union)
            CSV.write(outfile, combined)
        else
            println("  • $sample  ($nrows rows)  → writing $outfile")
            CSV.write(outfile, sub)
        end
    end
end

println("\n✅ Done – all files processed and split by sample")