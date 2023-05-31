VecF64 = Vector{Float64}

struct ControlPoint
    x::Float64
    w::Float64
end

struct Layer
    controlPoints::Vector{ControlPoint}
    h::Float64
end

struct Vessel
    layers::Vector{Layer}
end

function normalize_layer(layer::Layer)::Layer
    return Layer(sort(layer.controlPoints, by=cp -> cp.x), layer.h)
end

function normalize_vessel(vessel::Vessel)::Vessel
    return Vessel(sort(map(normalize_layer, vessel.layers), by=layer -> layer.h))
end

struct Interpolator
    pp::Any
    start_x::Float64
    start_w::Float64
    end_x::Float64
    end_w::Float64
end

function make_layer_interp(interp_fn, layer::Layer)::Interpolator
    @assert length(layer.controlPoints) >= 2

    xs = map(cp -> cp.x, layer.controlPoints)
    ws = map(cp -> cp.w, layer.controlPoints)

    return Interpolator(
        interp_fn(xs, ws),
        first(xs),
        first(ws),
        last(xs),
        last(ws)
    )
end

function call_interp(interp::Interpolator, x::Float64)::Union{Float64,Nothing}
    if (x < interp.start_x || x > interp.end_x)
        return nothing
    elseif x == interp.start_x
        return interp.start_w
    elseif x == interp.end_x
        return interp.end_w
    else
        rst = ppval(interp.pp, x)[1]
        if rst < 0.0
            return 0.0
        end

        return rst
    end
end

function call_interp_seq(interp_seq::Vector{Interpolator}, x::Float64)::Union{Float64,Nothing}
    for interp in interp_seq
        rst = call_interp(interp, x)
        if !isnothing(rst)
            return rst
        end
    end

    return nothing
end

function make_vessel_layers_interp(interp_fn, vessel::Vessel)::Vector{Tuple{Interpolator,Float64}}
    return map(vessel.layers) do layer
        return (make_layer_interp(interp_fn, layer), layer.h)
    end
end

function compute_layer_heights(vessel::Vessel)::VecF64
    return map(layer -> layer.h, vessel.layers)
end

function compute_xw(vessel::Vessel)::Tuple{Vector{VecF64},Vector{VecF64}}
    xs::Vector{VecF64} = map(layer -> map(cp -> cp.x, layer.controlPoints), vessel.layers)
    ws::Vector{VecF64} = map(layer -> map(cp -> cp.w, layer.controlPoints), vessel.layers)

    return (xs, ws)
end

function compute_forerear_xw(
    xs::Vector{VecF64},
    ws::Vector{VecF64}
)::Tuple{VecF64,VecF64,VecF64,VecF64}
    fore_xs = map(first, xs)
    rear_xs = map(last, xs)
    fore_ws = map(first, ws)
    rear_ws = map(last, ws)

    return (fore_xs, rear_xs, fore_ws, rear_ws)
end

function make_vessel_forerear_interp(
    interp,
    heights::VecF64,
    fore_xs::VecF64,
    rear_xs::VecF64,
    fore_ws::VecF64,
    rear_ws::VecF64
)::Tuple{Interpolator,Interpolator,Interpolator,Interpolator}
    min_h = first(heights)
    max_h = last(heights)

    return (
        Interpolator(interp(heights, fore_xs), min_h, first(fore_xs), max_h, last(fore_xs)),
        Interpolator(interp(heights, rear_xs), min_h, first(rear_xs), max_h, last(rear_xs)),
        Interpolator(interp(heights, fore_ws), min_h, first(fore_ws), max_h, last(fore_ws)),
        Interpolator(interp(heights, rear_ws), min_h, first(rear_ws), max_h, last(rear_ws))
    )
end

struct GenerationConfig
    section_length::Float64
    section_height::Float64
end

function compute_generated_heights(
    total_height::Float64,
    config::GenerationConfig
)::VecF64
    cnt = fld(total_height, config.section_height)
    cnt_h = cnt * config.section_height
    rst::VecF64 = collect(0.0:config.section_height:cnt_h)

    if cnt_h != total_height
        push!(rst, total_height)
    end

    return rst
end

function compute_hard_breaks(xs::Vector{VecF64})::VecF64
    return unique(sort(collect(Iterators.flatten(xs))))
end

function make_hard_break_interp(
    interp_fn::Any,
    vessel_layers_interp::Vector{Tuple{Interpolator,Float64}},
    hard_breaks::VecF64
)::Vector{Vector{Interpolator}}
    rst::Vector{Vector{Interpolator}} = []

    for hard_break in hard_breaks
        min_h::Union{Float64,Nothing} = nothing
        max_h::Union{Float64,Nothing} = nothing
        ws::Union{VecF64,Nothing} = nothing
        hs::Union{VecF64,Nothing} = nothing

        brk_interp::Vector{Interpolator} = []

        for (layer_interp, h) in vessel_layers_interp
            w = call_interp(layer_interp, hard_break)

            if isnothing(w)
                if !isnothing(min_h)
                    if length(ws) >= 2
                        interp_sect = Interpolator(
                            interp_fn(ws::VecF64, hs::VecF64),
                            min_h::Float64,
                            first(ws::VecF64),
                            max_h::Float64,
                            last(ws::VecF64)
                        )
                        push!(brk_interp, interp_sect)
                    end
                    min_h = nothing
                    max_h = nothing
                    ws = nothing
                    hs = nothing
                end
            else
                if isnothing(min_h)
                    min_h = h
                    ws = []
                    hs = []
                end
                max_h = h
                push!(ws, w)
                push!(hs, h)
            end
        end

        if !isnothing(min_h)
            if length(ws) >= 2
                interp_sect = Interpolator(
                    interp_fn(hs::VecF64, ws::VecF64),
                    min_h::Float64,
                    first(ws::VecF64),
                    max_h::Float64,
                    last(ws::VecF64)
                )
                push!(brk_interp, interp_sect)
            end
            min_h = nothing
            max_h = nothing
            ws = nothing
            hs = nothing
        end

        push!(rst, brk_interp)
    end

    return rst
end

function make_interp(
    interp_fn::Any,
    generated_heights::VecF64,
    interp_fore_x::Interpolator,
    interp_rear_x::Interpolator,
    interp_fore_w::Interpolator,
    interp_rear_w::Interpolator,
    hard_breaks::VecF64,
    hard_break_interp::Vector{Vector{Interpolator}}
)::Vector{Interpolator}
    rst::Vector{Interpolator} = []

    for h in generated_heights
        start_x = call_interp(interp_fore_x, h)
        end_x = call_interp(interp_rear_x, h)
        start_w = call_interp(interp_fore_w, h)
        end_w = call_interp(interp_rear_w, h)

        @assert !isnothing(start_x) && !isnothing(end_x)
        @assert !isnothing(start_w) && !isnothing(end_w)

        xs::VecF64 = [start_x]
        ws::VecF64 = [start_w]

        current_hardbreak::Int = 1
        while current_hardbreak <= length(hard_breaks) && hard_breaks[current_hardbreak] <= start_x
            current_hardbreak += 1
        end

        while current_hardbreak <= length(hard_breaks) && hard_breaks[current_hardbreak] <= end_x
            w = call_interp_seq(hard_break_interp[current_hardbreak], h)
            if !isnothing(w)
                push!(xs, hard_breaks[current_hardbreak])
                push!(ws, w::Float64)
            end
            current_hardbreak += 1
        end

        if last(xs) != end_x
            push!(xs, end_x)
            push!(ws, end_w)
        end

        push!(rst, Interpolator(
            interp_fn(xs, ws),
            start_x,
            start_w,
            end_x,
            end_w
        ))
    end

    return rst
end

struct DesignedVessel
    heights::VecF64
    interp::Vector{Interpolator}
    fore::VecF64
    rear::VecF64
end

function principalized_tessellate(
    interp_fn::Any,
    vessel::Vessel,
    cfg::GenerationConfig;
    plotName::Union{String,Nothing}=nothing,
    color::String="blue"
)::DesignedVessel
    vessel = normalize_vessel(vessel)

    layer_heights = compute_layer_heights(vessel)
    total_height = last(layer_heights)

    xs, ws = compute_xw(vessel)

    fore_xs, rear_xs, fore_ws, rear_ws = compute_forerear_xw(xs, ws)
    interp_fore_x, interp_rear_x, interp_fore_w, interp_rear_w = make_vessel_forerear_interp(
        interp_fn,
        layer_heights,
        fore_xs,
        rear_xs,
        fore_ws,
        rear_ws
    )

    generated_heights = compute_generated_heights(total_height, cfg)
    hard_breaks = compute_hard_breaks(xs)

    vessel_layer_interps = make_vessel_layers_interp(interp_fn, vessel)
    hard_break_interps = make_hard_break_interp(
        interp_fn,
        vessel_layer_interps,
        hard_breaks
    )

    interps = make_interp(
        interp_fn,
        generated_heights,
        interp_fore_x,
        interp_rear_x,
        interp_fore_w,
        interp_rear_w,
        hard_breaks,
        hard_break_interps
    )

    if !isnothing(plot)
        args = []
        for i in eachindex(interps)
            interp = interps[i]
            h = generated_heights[i]
            y1 = x -> call_interp(interp, x) / 2
            y2 = x -> -y1(x)

            xs1 = collect(interp.start_x:0.1:interp.end_x)
            ys1 = y1.(xs1)
            ys2 = y2.(xs1)

            if first(ys1) != 0.0
                pushfirst!(xs1, first(xs1) - 0.001)
                pushfirst!(ys1, 0.0)
                pushfirst!(ys2, 0.0)
            end

            if last(ys1) != 0.0
                push!(xs1, last(xs1) + 0.001)
                push!(ys1, 0.0)
                push!(ys2, 0.0)
            end

            push!(args, xs1, ys1, h, xs1, ys2, h)
        end

        figure(plotName, figsize=(8, 8))
        plot3(args...; color=color)
        daspect([1 1 1])
    end

    return DesignedVessel(
        generated_heights,
        interps,
        map(h -> call_interp(interp_fore_x, h), generated_heights),
        map(h -> call_interp(interp_rear_x, h), generated_heights)
    )
end
