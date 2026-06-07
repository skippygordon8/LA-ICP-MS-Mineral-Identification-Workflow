cd(@__DIR__)

using JSON3, StatGeochemBase, Isoplot, Measurements
using Plots

#Load standards (only once)

standardfilenames = ["NISTSRM616", "NISTSRM614", "NISTSRM612", "NISTSRM610"]
standardnames     = ["NIST616",    "NIST614",    "NIST612",    "NIST610"]

standards = Dict{String,Dict{String,Float64}}()

for i in eachindex(standardfilenames)
    d = Dict{String,Float64}()
    jobj = JSON3.read(read("$(standardfilenames[i]).json", String))
    for e in jobj.element
        d[e.label] = e.concentration_ppm
        d[e.label * "_sigma"] = e.concentration_ppm_stdev > 0 ?
                                e.concentration_ppm_stdev :
                                0.1 * e.concentration_ppm
    end
    standards[standardnames[i]] = d
end


colmap(data) = Dict(lowercase(String(k)) => k for k in keys(data))

# Build Calibrations (only once)

println("Building calibrations from standards file")   

sds    = importdataset("standards_file.csv", importas=:Tuple)    #Change file name to compiled standards dataset.
sdsmap = colmap(sds)

spot         = sdsmap["spot_id"]
standardrows = contains.(sds[spot], "NIST")
standardIDs  = sds[spot][standardrows] .|> s -> s[3:end-2]

# all CPS columns in the standards file
cps_keys = filter(k -> occursin("_cps", k) && !occursin("error", k),
                  collect(keys(sdsmap)))

# Dict:  element-base (e.g. "li7") =>  yorkfit result
fits = Dict{String,Any}()

for k in cps_keys
    base = replace(k, "_cps" => "")

    cps_sym     = sdsmap[k]
    cps_err_sym = sdsmap[base * "_cps_error"]

    element = replace(base, r"[0-9]+" => "") |> uppercasefirst

    y       = standardIDs .|> x -> standards[x][element]
    y_sigma = standardIDs .|> x -> standards[x][element * "_sigma"]

    x       = sds[cps_sym][standardrows]
    x_sigma = sds[cps_err_sym][standardrows]

    push!(y, 0);  push!(y_sigma, 1)
    push!(x, 0);  push!(x_sigma, 1)

 
    # fit + plot (once per element)

    yf = yorkfit(x, x_sigma, y, y_sigma)
    fits[base] = yf

    h = plot(title=element, xlabel="CPS", ylabel="PPM", framestyle=:box)
    scatter!(h, x, y, xerror=x_sigma, yerror=y_sigma)
    plot!(h, yf)

    plotname = "calibration_$(element).pdf"
    savefig(h, plotname)
    println("  • $element → $plotname")
end

println("Built $(length(fits)) calibrations.\n")


# APPLY CALIBRATIONS TO EACH DATA FILE


function process_all_files(fits)

    files = readdir(@__DIR__)

    datafiles = filter(f ->
        occursin("JGordon_DATA", f) && endswith(f, ".csv"), #change to match data file naming scheme
        files
    )

    println("Found $(length(datafiles)) valid data files.")

    for file in datafiles
        println("\nProcessing $file")

        ds    = importdataset(file, importas=:Tuple)
        dsmap = colmap(ds)

        ds_spot = dsmap["spot_id"]

        for (base, yf) in fits
            cps_key     = base * "_cps"
            cps_err_key = base * "_cps_error"
            ppm_key     = base * "_ppm"
            ppm_err_key = base * "_ppm_error"

            if !(haskey(dsmap, cps_key)     &&
                 haskey(dsmap, cps_err_key) &&
                 haskey(dsmap, ppm_key)     &&
                 haskey(dsmap, ppm_err_key))
                @warn "Skipping $base in $file (missing column)"
                continue
            end

            cps     = dsmap[cps_key]
            cps_err = dsmap[cps_err_key]
            ppm     = dsmap[ppm_key]
            ppm_err = dsmap[ppm_err_key]

            calc = Isoplot.line.(yf, ds[cps] .± ds[cps_err])

            ds[ppm]     .= value.(calc)
            ds[ppm_err] .= stdev.(calc)
        end

        # export

        ppm_keys = filter(k -> occursin("_ppm", k) && !occursin("error", k),
                          collect(keys(dsmap)))
        sort!(ppm_keys)

        export_cols = [ds_spot, [dsmap[k] for k in ppm_keys]...]

        prefix  = split(file, "_JGordon")[1]      # Change to match file naming scheme
        outfile = "$(prefix)_processed.csv"       

        exportdataset(ds, export_cols, outfile)

        println("✅ Saved dataset → $outfile")
    end
end

process_all_files(fits)