"""
# SeisPlot

Create plots of seismic data using [`Seis`](https://github.com/anowacki/Seis.jl)
and [`Plots`](https://github.com/JuliaPlots/Plots.jl)

## Plot types

### `plot`
Calling `plot` with either a single `Seis.Trace` or array of traces will produce
a set of plots, one for each trace.

### `section`
This produces a record section of an array of traces.
"""
module SeisPlot

using RecipesBase
using Seis

"Default maximum number of samples to plot at one time"
const MAX_SAMPLES = 30_000
"Default decimation behaviour.  May be useful to change if using a Plots backend which
does its own decimation or doesn't suffer performance issues."
const DECIMATE = Ref(true)

"""
    plot(traces::AbstractArray{<:Seis.Trace}; kwargs...) -> ::Plots.Plot

Plot a set of traces as a set of separate wiggles.

Additional options provided:

- `ylims`=`:all`: All traces have same amplitude limits
- `picks`: If `false`, do not plot picks.
"""
plot

# Recipe defining the above
@recipe function f(t::AbstractArray{<:Trace}; picks=true)
    ntraces = length(t)
    layout --> (ntraces, 1)
    # Plotting defaults
    legend --> false
    framestyle --> :box
    grid --> false
    linecolor --> :black
    linewidth --> 1
    # Set amplitude limits
    if get(plotattributes, :ylims, nothing) == :all
        plotattributes[:ylims] = extrema(vcat(trace.(t)...))
    end
    # Time limits
    if !haskey(plotattributes, :xlims)
        xlims = (minimum(first.(times.(t))), maximum(last.(times.(t))))
    else
        xlims = plotattributes[:xlims]
    end
    # Plot traces
    for i in eachindex(t)
        @series begin
            subplot := i
            xticks --> (i < ntraces ? nothing : :auto)
            times(t[i]), trace(t[i])
        end
    end

    # Picks
    if picks
        # Pick lines
        seriestype := :vline
        linecolor := :blue
        linewidth --> 1
        annot_params = (10, :blue, :left, :bottom)
        for i in eachindex(t)
            @series begin
                xticks --> (i < ntraces ? nothing : :auto)
                guidefontcolor := :blue
                subplot := i
                [p.time for p in Seis.picks(t[i]) if xlims[1] <= p.time <= xlims[2]]
            end
        end

        # Pick labels: workaround (see https://discourse.julialang.org/t/annotations-and-line-widths-in-plots-jl-heatmaps/4259/3)
        seriestype := :scatter
        markerstrokealpha := 0.0
        # seriesalpha := 0.0
        markeralpha := 0.0
        primary := false
        annotation_params = (8, :left, :bottom, :blue)
        for i in eachindex(t)
            @series begin
                xticks --> (i < ntraces ? nothing : :auto)
                subplot := i
                # FIXME: Update to a better way of implementing this: Main.Plots may break
                series_annotations := [Main.Plots.text.(coalesce(p.name, ""), annotation_params...)
                                       for p in Seis.picks(t[i])
                                       if xlims[1] <= p.time <= xlims[2]]
                x = [p.time for p in Seis.picks(t[i]) if xlims[1] <= p.time <= xlims[2]]
                y = get(plotattributes, :ylims, extrema(trace(t[i])))[1]
                x, repeat([y], length(x))
            end
        end
    end
end

"""
    section(traces, y_values=distance_deg.(traces); kwargs...) -> ::Plots.plot

Plot a record section for a collection of `Trace`s.  By default, the y-axis is
epicentral distance in degrees, but specify one of the following for `y_values`
to change this:

- A `Function` which will be applied to the traces, e.g., `distance_km`.
- A `Symbol` which corresponds to a field of the `Trace`s' `meta` dictionary.
- A set of values (i.e., `AbstractArray`)
- One of the following strings:
  - "index" for trace index (`i:length(traces)`)

Additional options provided via keyword arguments:

- `align`: Set to a `String` to align on the first pick with this name.
           Set to an array of values to align on the value for each trace.
- `decimate`: If `false`, do not perform downsampling of traces for plotting.
           Defaults to `true`.
- `max_samples`: Control the maximum number of samples to display at one time
           in order to make plotting quicker.  Set `decimate` to `false` to turn
           this off.
"""
section

@userplot Section

@recipe function f(sec::Section; align=nothing, decimate=DECIMATE[], max_samples=MAX_SAMPLES)
    # Arguments
    t = sec.args[1]
    y_values = if length(sec.args) >= 2
        sec.args[2]
    else
        ylabel --> "Distance / °"
        distance_deg.(t)
    end
    t isa AbstractArray{<:Trace} ||
        throw(ArgumentError("section requires an array of `Seis.Trace`s"))
    decimate isa Bool || throw(ArgumentError("decimate must be `true` or `false`"))
    # Defaults
    linecolor --> :black
    linewidth --> 1
    legend --> false
    framestyle --> :box
    grid --> false
    xlabel --> "Time / s"
    # Time shifts
    shifts = if align == nothing
        zeros(length(t))
    elseif align isa String || align isa Regex
        shifts = try
            [first(picks(tt, align)).time for tt in t]
        catch err
            error("Error finding pick '$align' for every trace.  (Error: $err)")
        end
        shifts
    elseif align isa AbstractVector
        length(align) == length(t) ||
            throw(ArgumentError("`align` must have same length as number of traces"))
        align
    else
        throw(ArgumentError("`align` must be a pick name or vector of shifts"))
    end
    # Y-axis
    y_shifts = if y_values isa Function
        y_values.(t)
    elseif y_values isa Symbol
        getproperty.(t.meta, y_values)
    elseif y_values isa AbstractArray
        length(y_values) == length(t) ||
            throw(ArgumentError("Length of y values ($(length(y_values))) does not " *
                                "match number of traces ($(length(t)))"))
        y_values
    elseif y_values isa AbstractString
        # Custom
        if y_values == "index"
            1:length(t)
        else
            throw(ArgumentError("Unrecognised y axis name '$y_values'"))
        end
    end
    # Get decimation value
    t1, t2 = get(plotattributes, :xlims, (-Inf, Inf))
    ndecimate = decimate ? decimation_value(t, shifts, t1, t2, max_samples) : 1
    # Traces
    time = [times(tt)[1:ndecimate:end] .- shift for (tt,shift) in zip(t, shifts)]
    maxval = maximum([maximum(abs.(tt)) for tt in trace.(t)])
    traces = [(trace(tt)./maxval .+ y)[1:ndecimate:end] for (tt, y) in zip(t, y_shifts)]
    # Plot
    @series begin
        time, traces
    end
end

"""
    decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples) -> n
    
Return the decimation value `n` which ensures that no more than `max_samples`
are contained within the time window `t1`-`t2` s relative to the values in `shifts`.
"""
function decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples)
    n = 1
    np = sum(num_points.(t, shifts .+ t1, shifts .+ t2))
    while np÷n > max_samples
        n += 1
    end
    n
end

"""
    num_points(t::Trace, t1, t2) -> np

Return the number of points `np` between `t1` and `t2` seconds relative to the
origin time of the trace `t`.
"""
num_points(t::Trace, t1, t2) = sum(t1 .<= times(t) .<= t2)

end # module
