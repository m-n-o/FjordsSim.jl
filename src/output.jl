using CairoMakie: Axis, Figure, Colorbar, Observable, Reverse, record, heatmap!, @lift, lines!, axislegend
using FileIO: save
using Oceananigans.Fields: Field, interior, compute!

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_1d_phys(T, S, z, times, folder, x, y)
    fig = Figure(size = (1000, 400), fontsize = 20)

    axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", xticks = (0:30:times[end]/days), xtickformat = "{:.0f}")

    Axis(fig[1, 1]; title = "T, ⁰C", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, x, y, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[1, 2], hmT)

    Axis(fig[2, 1]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, x, y, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 2], hmS)

    save(joinpath(folder, "1d_phys.png"), fig)
end

map_axis_kwargs = (xlabel = "Grid points, eastward direction", ylabel = "Grid points, northward direction")
transect_axis_kwargs = (xlabel = "Grid points, eastward direction", ylabel = "z (m)")
framerate = 12

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_surface_speed(u, v, Nz, times, folder; colorrange = (0, 0.5), colormap = :deep)
    Nt = length(times)
    iter = Observable(Nt)

    ## Speed
    si = @lift begin
        s = Field(sqrt(u[$iter]^2 + v[$iter]^2))
        compute!(s)
        s = interior(s, :, :, Nz)
        s[s.==0] .= NaN
        s
    end

    fig = Figure(size = (1000, 400))

    title = @lift "Surface speed at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, si, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface speed (ms⁻¹)")
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "surface_speed.gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_bottom_tracer(
    variable,
    var_name,
    bottom_z,
    times,
    folder;
    colorrange = (-1, 350),
    colormap = :turbo,
    figsize = (1000, 400),
)

    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        Ti = [variable[i, j, bottom_z[i, j], $iter] for i = 1:size(variable, 1), j = 1:size(variable, 2)]
        Ti[Ti.==0] .= NaN
        Ti
    end

    title = @lift "bottom $(var_name), mmol/m³ at " * prettytime(times[$iter])
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name), mmol/m³")

    Nt = length(times)
    record(fig, joinpath(folder, "bottom_$(var_name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_horizontal_tracer(tracer, times, folder, name, label; colorrange = (-1, 30), colormap = :magma, iz = 10)
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        if typeof(tracer) == Array{Float64, 4}
            Ti = tracer[:,:,iz,$iter]
        else
            Ti = interior(tracer[$iter], :, :, iz)
        end
        Ti[Ti.==0] .= NaN
        Ti
    end

    title = @lift label * " at " * prettytime(times[$iter])
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_vertical_tracer(tracer, depth, iy, times, folder, name, label; colorrange = (-1, 30), colormap = :magma)

    xs = 1:size(tracer)[1] # get x-values for x-axis
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        if typeof(tracer) == Array{Float64, 4}
            Ti = tracer[:,iy,:,$iter]
        else
            Ti = interior(tracer[$iter], :, iy, :)
        end
        Ti[Ti.==0] .= NaN
        Ti
    end

    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_vertical_tracer_points(
    tracer,
    depth,
    indices::Vector{Tuple{Int, Int}},  # List of (ix, iy) pairs
    times,
    folder,
    name,
    label;
    colorrange = (-1, 30),
    colormap = :magma,
)

    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        # Collect slices from all given (ix, iy) points
        slices = [interior(tracer[$iter], ix, iy, :) for (ix, iy) in indices]
        Ti = hcat(slices...)  # Stack slices
        Ti[Ti.==0] .= NaN
        Ti = Ti'
        Ti
    end

    xs = 1:size(indices)[1] # get x-values for x-axis
    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function record_vertical_diff(
    tracer,
    depth,
    iy,
    times,
    folder,
    name,
    label;
    colorrange = (-1, 30),
    colormap = :magma,
)

    xs = 1:size(tracer)[1] # get x-values for x-axis
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        Ti = tracer[:, iy, :, $iter]
        Ti[Ti.==0] .= NaN
        Ti
    end

    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_ztime(NUT, O₂, O₂_relative, PHY, HET, T, DOM, POM, S, i, j, times, z, folder)

    fig = Figure(size = (1500, 1000), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:50:times[end]),
        xtickformat = "{:.0f}", #   values -> ["$(value)kg" for value in values]     
    )

    axNUT = Axis(fig[1, 1]; title = "NUT, mmolN/m³", axis_kwargs...)
    hmNUT = heatmap!(times / days, z, interior(NUT, i, j, :, :)', colormap = Reverse(:cherry))
    Colorbar(fig[1, 2], hmNUT)

    axOXY = Axis(fig[1, 3]; title = "O₂, mmol/m³", axis_kwargs...)
    hmOXY = heatmap!(times / days, z, interior(O₂, i, j, :, :)', colormap = :turbo)
    Colorbar(fig[1, 4], hmOXY)
    
    axOXY_rel = Axis(fig[1, 5]; title = "O₂ saturation, %", axis_kwargs...)
    hmOXY_rel = heatmap!(times / days, z, O₂_relative[i, j, :, :]', colormap = :gist_stern)
    Colorbar(fig[1, 6], hmOXY_rel)

    axPHY = Axis(fig[2, 1]; title = "PHY, mmolN/m³", axis_kwargs...)
    hmPHY = heatmap!(times / days, z, interior(PHY, i, j, :, :)', colormap = Reverse(:cubehelix)) #(:davos10))
    Colorbar(fig[2, 2], hmPHY)

    axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
    hmHET = heatmap!(times / days, z, interior(HET, i, j, :, :)', colormap = Reverse(:afmhot))
    Colorbar(fig[2, 4], hmHET)

    axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, i, j, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 6], hmT)

    axDOM = Axis(fig[3, 1]; title = "DOM, mmolN/m³", axis_kwargs...)
    hmDOM = heatmap!(times / days, z, interior(DOM, i, j, :, :)', colormap = Reverse(:CMRmap)) #(:devon10))
    Colorbar(fig[3, 2], hmDOM)

    axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
    hmPOM = heatmap!(times / days, z, interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
    Colorbar(fig[3, 4], hmPOM)

    axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, i, j, :, :)', colormap = :viridis)
    Colorbar(fig[3, 6], hmS)

    save(joinpath(folder, "ztime_$(i)_$(j).png"), fig)
   
    @info "VARIABLES ztime_$(i)_$(j) plots made"
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_bottom_tracer(tracer, name, bottom_z, times, plot_day, folder; colorrange = (0, 0.5), colormap = :turbo, )

    #bottom_tracer = [tracer[i, j, bottom_z[i, j], time] for i = 1:size(tracer, 1), j = 1:size(tracer, 2)]
    #fig = Figure(size = (1000, 400), fontsize = 20)

    nday = plot_day*round(Int, (length(times)/365))
    #   nday = plot_day*round(Int, ( times[end]/365))
    z_tracer = [tracer[i, j, bottom_z[i, j], nday] == 0 ? NaN : tracer[i, j, bottom_z[i, j], nday] for i in 1:size(tracer, 1), j in 1:size(tracer, 2)]
   # z_tracer = [tracer[i, j, bottom_z[i, j], nday] for i = 1:size(tracer, 1), j = 1:size(tracer, 2)]
   #z_tracer = [tracer[i, j, bottom_z[i, j], plot_day] for i = 1:size(tracer, 1), j = 1:size(tracer, 2)]    
    fig = Figure(size = (1000, 500), fontsize = 20)
    axis_kwargs = (xlabel = "Grid points, East ", ylabel = "Grid points, North")
    axTRAC = Axis(fig[2, 1]; title = "$(name), mmol/m³, day $(nday) ", axis_kwargs...) 
    hmTRAS = heatmap!([i for i = 1:size(tracer, 1)], [j for j = 1:size(tracer, 2)], z_tracer, colorrange = colorrange, colormap = colormap, nan_color = :gray)
    Colorbar(fig[2, 2], hmTRAS)
    save(joinpath(folder, "bottom_$(name)_day_$(plot_day).png"), fig)   
    @info "bottom_$(name)_day_$(plot_day) plot made"    
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_surface_tracer(tracer, name, iz, times, plot_day, folder; colorrange = (0, 0.5), colormap = :turbo, )

   nday = plot_day*round(Int, (length(times)/365))
#   nday = plot_day*round(Int, ( times[end]/365))
   z_tracer = [tracer[i, j, iz, nday] == 0 ? NaN : tracer[i, j, iz, nday] for i in 1:size(tracer, 1), j in 1:size(tracer, 2)]

   fig = Figure(size = (1000, 500), fontsize = 20)
    axis_kwargs = (xlabel = "Grid points, East ", ylabel = "Grid points, North")
    axTRAC = Axis(fig[2, 1]; title = "$(name), mmol/m³, day $(nday) ", axis_kwargs...)
    hmTRAS = heatmap!([i for i = 1:size(tracer, 1)], [j for j = 1:size(tracer, 2)], z_tracer, colorrange = colorrange, colormap = colormap, nan_color = :grey) #:gray)
    Colorbar(fig[2, 2], hmTRAS)

    save(joinpath(folder, "surface_$(name)_day_$(plot_day).png"), fig)
    @info "surface_$(name)_day_$(plot_day) plot made"    
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_ratio_above_thresh(variable, name, bottom_z, left_x, right_x, times, folder)

    # bottom_z evaluation
    variable = variable[left_x:right_x,:,:,:]
    #@info size(variable)

    bottom_z = bottom_z[left_x:right_x,:]
    #@info size(bottom_z)

    # bottom variable Matrix evaluation
    variable_bot = Array{eltype(variable)}(undef, size(variable, 1), size(variable, 2), size(variable, 4))
    for i in 1:size(variable, 1)
        for j in 1:size(variable, 2)
            variable_bot[i, j, :] = variable[i, j, bottom_z[i, j], :]
        end
    end

    land_cells = length(variable_bot[:,:,1][variable_bot[:,:,1] .< 10])
    water_cells = length(bottom_z) - land_cells
    
    under_thresh(thresh) = [length(variable_bot[:,:,i][variable_bot[:,:,i] .< thresh]) - land_cells for i in 1:size(variable_bot, 3)]
    
    fig = Figure(size = (400, 400), fontsize = 20) #    fig = Figure(size = (1000, 400), fontsize = 20)
    
     axis_kwargs = (
    xlabel = "Time (days)",
    ylabel = "% of bottom area",
    xticks = (0:60:times[end]/days),
    xtickformat = "{:.0f}",
    limits = ((nothing, nothing), (0, 50)),   # y-axis from 0 to 50
    #aspect = DataAspect()                    # equal scaling for both axes
)
    Axis(fig[1, 1]; axis_kwargs...)
    
    lines!(times / days, 100*under_thresh(90)/water_cells, colormap = Reverse(:RdYlBu), label="$(name)< 90 μM")
    lines!(times / days, 100*under_thresh(15)/water_cells, colormap = Reverse(:RdYlBu), label="$(name)< 15 μM")
    lines!(times / days, 100*under_thresh(1)/water_cells, colormap = Reverse(:RdYlBu), label="$(name)< 1 μM")
    axislegend()

    save(joinpath(folder, "ratio_$(name)_above_thresh_$(left_x)-$(right_x).png"), fig)
    @info "ratio_$(name)_above_thresh_$(left_x)-$(right_x) plot made"
end