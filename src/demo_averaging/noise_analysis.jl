cv(x) = std(x) / mean(x)
fano(x) = var(x) / mean(x)
function noiseplots(;n_upsample=4)
    r = LinRange(0,1,201)
    seam_cell_positions_over_time = map(eachindex(r)) do j
        _model = avg_models[j]
        seam_cell_pts(_model, n_upsample)
    end
    seam_cell_positions_over_time_matrix = stack(seam_cell_positions_over_time) .* voxel_size
    cfc = CoordinateTransformations.CylindricalFromCartesian()
    r = cfc.(seam_cell_positions_over_time_matrix) .|> (x->x.r)
    z = cfc.(seam_cell_positions_over_time_matrix) .|> (x->x.z)
    fig = Figure()
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])
    ax3 = Axis(fig[1,2])
    ax4 = Axis(fig[2,2])
    ax5 = Axis(fig[1,3])
    ax6 = Axis(fig[2,3])

    plot!(ax1, (var.(eachrow(r)))[1:end÷2])
    ax1.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax1.xlabel = "Seam Cells"
    ax1.ylabel = "Variance (σ^2)"
    ax1.title = "Variance of R for Seam Cells"

    plot!(ax2, (var.(eachrow(z)))[1:end÷2])
    ax2.xlabel = "Seam Cells"
    ax2.ylabel = "Variance (σ^2)"
    ax2.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax2.title = "Variance of Z for Seam Cells"

    plot!(ax3, (cv.(eachrow(r)))[1:end÷2])
    ax3.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax3.xlabel = "Seam Cells"
    ax3.ylabel = "Coefficient of Variation (σ/μ)"
    ax3.title = "CV of R for Seam Cells"

    plot!(ax4, (cv.(eachrow(z)))[1:end÷2])
    ax4.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax4.xlabel = "Seam Cells"
    ax4.ylabel = "Coefficient of Variation (σ/μ)"
    ax4.title = "CV of Z for Seam Cells"

    plot!(ax5, (fano.(eachrow(r)))[1:end÷2])
    ax5.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax5.xlabel = "Seam Cells"
    ax5.ylabel = "Fano Factor (σ^2/μ)"
    ax5.title = "Fano Factor of R for Seam Cells"

    plot!(ax6, (fano.(eachrow(z)))[1:end÷2])
    ax6.xticks = (1:11, String.(avg_models[1].names[1:2:end]))
    ax6.xlabel = "Seam Cells"
    ax6.ylabel = "Fano Factor (σ^2/μ)"
    ax6.title = "Fano Factor of Z for Seam Cells"


    return fig
end