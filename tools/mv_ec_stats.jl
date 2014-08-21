using JSON
using Glob

raw_ec_records = [JSON.parsefile(a) for a in ARGS]

# 3x3 arrays of counts, 2 symbols, 4 levels
ec_stats = zeros(9, 2, 4)

for frames in raw_ec_records
    for frame in frames
        key_levels = ["mvf-l1" => 1, "mvf-l2" => 2, "mvf-l3" => 3, "mvf-l4" => 4]
        for key in keys(key_levels)
            if haskey(frame, key)
                for r in frame[key]
                    sym = r[1] + 1
                    ctx = r[3] + 1
                    ec_stats[ctx, sym, key_levels[key]] += 1
                end
            end
        end
    end
end

for level in 1:4
    totals = ec_stats[:,1,level] .+ ec_stats[:,2,level]
    totals[find(isnan,totals)] = 0.0
    probs = ec_stats[:,1,level] ./ totals
    probs[find(isnan,probs)] = 0.0
    probs = int(round(probs * 32768))
    println("LEVEL $level:")
    println("Probabilties:")
    println(probs)
    println("Totals:")
    println(totals)
    println()
end
