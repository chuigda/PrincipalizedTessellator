include("./vessel.jl")

function test()
    v::Vessel = Vessel([
        Layer([
            ControlPoint(0, 0),
            ControlPoint(1, 1),
            ControlPoint(2, 2),
            ControlPoint(4, 3),
            ControlPoint(8, 4),
            ControlPoint(16, 5),
            ControlPoint(32, 6)
        ], 4),
        Layer([
            ControlPoint(1 + 4, 0),
            ControlPoint(2 + 4, 1),
            ControlPoint(4 + 4, 2),
            ControlPoint(8 + 4, 3),
            ControlPoint(16 + 4, 4),
            ControlPoint(31, 5)
        ], 0)
    ])

    cfg::GenerationConfig = GenerationConfig(3.0, 0.5)

    principalized_tessellate(spline, v, cfg, plotName="spline", color="blue")
    principalized_tessellate(pchip, v, cfg, plotName="pchip", color="red")

    return nothing
end

test()
