import JSON

include("./vessel.jl")

function layer_from_json_object(json::Dict{String,Any})::Layer
    return Layer(
        map(json["controlPoints"]) do cp
            ControlPoint(
                parse(Float64, cp["x"]),
                parse(Float64, cp["w"])
            )
        end,
        parse(Float64, json["h"])
    )
end

function test()
    open("North Darkota.json") do f
        cfg::GenerationConfig = GenerationConfig(3.0, 1.0)

        s = read(f, String)
        json = JSON.parse(s)

        v::Vessel = Vessel([
            layer_from_json_object(json["bottom"])
            map(layer_from_json_object, json["layers"])
            layer_from_json_object(json["deck"])
        ])

        principalized_tessellate(spline, v, cfg, plotName="spline", color="blue")
        principalized_tessellate(pchip, v, cfg, plotName="pchip", color="red")
    end
end

test()
