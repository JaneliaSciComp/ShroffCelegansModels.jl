function make_angle_heatmap(smodels)
    angles = map(smodels) do smodel
        if !ismissing(smodel)
            θ = stretched_analysis(smodel)
            θ .- circ_mean(θ).μ
        else
            missing
        end
    end
    valid_angles = skipmissing(angles) |> collect
    angle_heatmap = valid_angles .|>
        (x->interpolate(LinRange(0,1,length(x)), x, BSplineOrder(4))) .|>
        (x->x.(r)) |> x->hcat(x...) .|>
        (x->mod(x, 2π)) |>
        (x->heatmap(x, colormap=:cyclic_mygbm_30_95_c78_n256))
    Colorbar(angle_heatmap.figure[:, end+1], angle_heatmap.plot)
    angle_heatmap.figure.content[2].label = "Angle (radians)"
    angle_heatmap.axis.ylabel = "Timepoints"
    angle_heatmap.axis.xlabel = "A-P Distance"
    
    angle_matrix = valid_angles .|> (x->interpolate(LinRange(0,1,length(x)), x, BSplineOrder(4))) .|> (x->x.(r)) |> x->hcat(x...) .|> (x->mod(x, 2π))

    scatter!(Axis(angle_heatmap.figure[1,3]), circ_var.(eachcol(angle_matrix)) .|> (x-> x.S), 1:size(angle_matrix,2))
    scatter!(Axis(angle_heatmap.figure[2,1]), circ_var.(eachrow(angle_matrix)) .|> (x-> x.S))

    return angle_heatmap
end