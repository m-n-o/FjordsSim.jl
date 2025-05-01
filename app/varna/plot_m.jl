using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie
using FjordsSim:
    plot_1d_phys,
    extract_z_faces,
    record_vertical_tracer,
    record_surface_speed,
    record_horizontal_tracer,
    plot_ztime,
    record_bottom_tracer,
    record_vertical_diff,
    plot_bottom_tracer,
    plot_surface_tracer,
    plot_ratio_above_thresh,
    oxygen_saturation

# garg = jldopen("$filename.jld2")["grid"]["underlying_grid"]
# grid_from_file = LatitudeLongitudeGrid{Bounded,Bounded,Bounded}(CPU(),
#     garg["Nx"], garg["Ny"], garg["Nz"],
#     garg["Hx"], garg["Hy"], garg["Hz"],
#     garg["Lx"], garg["Ly"], garg["Lz"],
#     garg["Δλᶠᵃᵃ"], garg["Δλᶜᵃᵃ"], garg["λᶠᵃᵃ"], garg["λᶜᵃᵃ"],
#     garg["Δφᵃᶠᵃ"], garg["Δφᵃᶜᵃ"], garg["φᵃᶠᵃ"], garg["φᵃᶜᵃ"], garg["z"],
#     garg["Δxᶠᶜᵃ"], garg["Δxᶜᶠᵃ"], garg["Δxᶠᶠᵃ"], garg["Δxᶜᶜᵃ"],
#     garg["Δyᶠᶜᵃ"], garg["Δyᶜᶠᵃ"],
#     garg["Azᶠᶜᵃ"], garg["Azᶜᶠᵃ"], garg["Azᶠᶠᵃ"], garg["Azᶜᶜᵃ"], garg["radius"])

include("setup.jl");

sim_setup = setup_region_3d();
grid_args = merge(sim_setup.grid_args, (arch = CPU(),))
grid_from_setup = sim_setup.grid_callable(grid_args...).underlying_grid
z = znodes(grid_from_setup, Center())
Nz = grid_from_setup.Nz

folder = joinpath(homedir(), "FjordsSim_results", "varna")
#filename = joinpath(folder, "pn30d40n00d00")  #ok
#filename = joinpath(folder, "pn30d40n01d02")  #ok
filename = joinpath(folder, "pn60d80n01d02")  #ok
#filename = joinpath(folder, "n30d40n01n02")  #ok
#filename = joinpath(folder, "n30d40n00_5d01")  #ok
#filename = joinpath(folder, "n30d40n00d00")  #ok

T = FieldTimeSeries("$filename.jld2", "T")
S = FieldTimeSeries("$filename.jld2", "S")
u = FieldTimeSeries("$filename.jld2", "u")
v = FieldTimeSeries("$filename.jld2", "v")
O₂ = FieldTimeSeries("$filename.jld2", "O₂")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
PHY = FieldTimeSeries("$filename.jld2", "P")
HET = FieldTimeSeries("$filename.jld2", "HET")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
POM = FieldTimeSeries("$filename.jld2", "POM")
#C = FieldTimeSeries("$filename.jld2", "C")
times = T.times

    @info "BGH arrays loaded"

# Calculate additional fields

Pressure = similar(interior(T)[:,:,:,1])
for k in 1:20
    Pressure[:,:,k] .= 100 + 9.8*k
end

O₂_saturation = similar(T)

for n in 1:length(T.times)
    T_slice = interior(T)[:, :, :, n]
    S_slice = interior(S)[:, :, :, n]
    P_slice = Pressure
    O₂_saturation[:, :, :, n] =@. oxygen_saturation(T_slice, S_slice, P_slice)
    compute!(O₂_saturation) 
end

O₂_relat = 100 .* O₂./O₂_saturation
O₂_relat[O₂_relat.==0] .= NaN

    @info "O2 % calculated"

bottom_z = ones(Int, size(O₂, 1), size(O₂, 2))
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = 1:size(O₂, 3)
            if O₂[i, j, k, 1] .!= 0
                bottom_z[i, j] = k
                break
                if k == 20
                    bottom_z[i, j] = 20
                end
            end
            if k == 20
                bottom_z[i, j] = 20
            end
        end
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  vertical distributions changes in a point
# ~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ztime(NUT, O₂, O₂_relat, PHY, HET, T, DOM, POM, S, 84, 14, times, z, folder)
plot_ztime(NUT, O₂, O₂_relat, PHY, HET, T, DOM, POM, S, 18, 18, times, z, folder)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# maps at  selected time 
# ~~~~~~~~~~~~~~~~~~~~~~~~~
plotday = 224
plot_surface_tracer(NUT, "NUT", Nz, times, plotday, folder; colorrange = (0, 40),  colormap = Reverse(:cherry), )
plot_surface_tracer(PHY, "PHY", Nz, times, plotday, folder; colorrange = (0, 7),   colormap = Reverse(:cubehelix), )
plot_surface_tracer(HET, "HET", Nz, times, plotday, folder; colorrange = (0, 5),  colormap = Reverse(:afmhot), )
plot_surface_tracer(POM, "POM", Nz, times, plotday, folder; colorrange = (0, 5),  colormap = Reverse(:greenbrownterrain), )
plot_surface_tracer(DOM, "DOM", Nz, times, plotday, folder; colorrange = (0, 40),  colormap = Reverse(:CMRmap), )
plot_surface_tracer(O₂,   "O₂", Nz, times, plotday, folder; colorrange = (0, 350), colormap = :turbo, )
plot_surface_tracer(O₂_relat, "O₂_%", Nz, times, plotday, folder; colorrange = (0, 120),colormap = :gist_stern, )

plot_bottom_tracer(NUT, "NUT", bottom_z, times, plotday, folder; colorrange = (0, 50),  colormap = Reverse(:cherry), )
plot_bottom_tracer(PHY, "PHY", bottom_z, times, plotday, folder; colorrange = (0, 7),   colormap = Reverse(:cubehelix), )
plot_bottom_tracer(HET, "HET", bottom_z, times, plotday, folder; colorrange = (0, 50),  colormap = Reverse(:afmhot), )
plot_bottom_tracer(POM, "POM", bottom_z, times, plotday, folder; colorrange = (0, 50),  colormap = Reverse(:greenbrownterrain), )
plot_bottom_tracer(DOM, "DOM", bottom_z, times, plotday, folder; colorrange = (0, 50),  colormap = Reverse(:CMRmap), )
plot_bottom_tracer(O₂,   "O₂", bottom_z, times, plotday, folder; colorrange = (0, 350), colormap = :turbo, )
plot_bottom_tracer(O₂_relat, "O₂_%", bottom_z, times, plotday, folder; colorrange = (0, 120),colormap = :gist_stern, )

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  Oxygen depletion plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# here we prescribe for the Varna case 3 intervals for "j":
# Beloslav: 1-28; Varna Lake: 42-120 and Varna Bay; 128-147
plot_ratio_above_thresh(O₂, "O₂", bottom_z, 1, 28, times, folder)
plot_ratio_above_thresh(O₂, "O₂", bottom_z, 42, 120, times, folder)
plot_ratio_above_thresh(O₂, "O₂", bottom_z, 128, 147, times, folder)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  Old and misc
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_variable(T, "temperature surface", Nz, T.times, folder, (1000, 400); colorrange = (-1, 20))
# record_variable(S, "salinity surface", Nz, S.times, folder, (1000, 400); colorrange = (0, 30))
# record_variable(u, "u velocity surface", Nz, u.times, folder, (1000, 400); colorrange = (-1, 1))
# record_variable(x_momentum, "x momentum", 1, x_momentum.times, folder, (1000, 400); colorrange = (-0.1, 0.1))
# record_variable(v, "v velocity surface", Nz, v.times, folder, (1000, 400); colorrange = (-1, 1))
# record_variable(y_momentum, "y momentum", 1, y_momentum.times, folder, (1000, 400); colorrange = (-0.1, 0.1))
# record_bottom_tracer(O₂, "Oxygen", Nz, O₂.times, folder)
# println(keys(grid["underlying_grid"]))
# println(grid["underlying_grid"]["Δyᶠᶜᵃ"])
 record_surface_speed(u, v, Nz, times, folder)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  HORIZONTAL animations bottom
# ~~~~~~~~~~~~~~~~~~~~~~~~~
record_bottom_tracer(
    O₂,
    "O₂",
    bottom_z,
    times,
    folder,
    colorrange = (0, 350),
    colormap = :turbo,
    figsize = (1000, 400)
)

record_bottom_tracer(
    HET,
    "HET",
    bottom_z,
    times,
    folder,
    colorrange = (0, 50),
    colormap =  Reverse(:afmhot),
    figsize = (1000, 400)
)

record_bottom_tracer(
    O₂_relat,
    "O₂_relat",
    bottom_z,
    times,
    folder,
    colorrange = (0, 110),
    colormap = :gist_stern,
    figsize = (1000, 400)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#  HORIZONTAL animations surface
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_horizontal_tracer(
#     C,
#     times,
#     folder,
#     "Contsurf",
#     "Contaminant (% of max. concentration)",
#     colorrange = (0, 100),
#     colormap = :matter,
#     iz = Nz,
# )

# record_horizontal_tracer(
#     T,
#     times,
#     folder,
#     "Tsurf",
#     "Temperature (°C)",
#     colorrange = (5, 40),
#     colormap = Reverse(:RdYlBu),
#     iz = Nz,
# )

# record_horizontal_tracer(
#     S,
#     times,
#     folder,
#     "Ssurf",
#     "Salinity (PSU)",
#     iz = Nz,
#     colorrange = (0, 17),
#     colormap = :viridis,
# )

 record_horizontal_tracer(
     O₂,
     times,
     folder,
     "O2surf",
     "Dissolved oxygen (μM)",
     colorrange = (0, 350),
     colormap = :turbo,
     iz = Nz,
 )

record_horizontal_tracer(
    O₂_relat,
    times,
    folder,
    "O₂_relat_surf",
    "%",
    iz = 20,
    colorrange = (0, 110),
    colormap = :gist_stern,
)

 record_horizontal_tracer(
     NUT,
     times,
     folder,
     "NUTsurf",
     "Nutrients (μM N)",
     colorrange = (0, 25),
     colormap = Reverse(:cherry),
     iz = Nz,
 )

 record_horizontal_tracer(
     PHY,
     times,
     folder,
     "PHYsurf",
     "Phytoplankton (μM N)",
     colorrange = (0, 5),
     colormap = Reverse(:cubehelix),
     iz = Nz,
 )

 record_horizontal_tracer(
     HET,
     times,
     folder,
     "HETsurf",
     "Heterotrophs (μM N)",
     colorrange = (0, 5),
     colormap = Reverse(:afmhot),
     iz = Nz,
 )

 record_horizontal_tracer(
     DOM,
     times,
     folder,
     "DOMsurf10",
     "DOM (μM N)",
     colorrange = (0, 50),
     colormap = Reverse(:CMRmap),
     iz = Nz,
 )

 record_horizontal_tracer(
     POM,
     times,
     folder,
     "POMsurf10",
     "POM (μM N)",
     colorrange = (0, 50),
     colormap = Reverse(:greenbrownterrain),
     iz = Nz,
 )


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# VERTICAL animations
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_vertical_tracer(
#     T,
#     z,
#     18,
#     times,
#     folder,
#     "Tprofile",
#     "Temperature (°C)",
#     colorrange = (5, 40),
#     colormap = Reverse(:RdYlBu),
# )

# record_vertical_tracer(
#     u,
#     z,
#     18,
#     times,
#     folder,
#     "Uprofile",
#     "u velocity component (ms⁻¹)",
#     colorrange = (-0.5, 0.5),
#     colormap = :deep,
# )

# record_vertical_tracer(S, z, 18, times, folder, "Sprofile", "Salinity (PSU)", colorrange = (0, 17), colormap = :viridis)

# record_vertical_tracer(
#     PHY,
#     z,
#     18,
#     times,
#     folder,
#     "PHYprofile",
#     "Phytoplankton (μM N)",
#     colorrange = (0, 5),
#     colormap = Reverse(:cubehelix),
# )

# record_vertical_tracer(
#     HET,
#     z,
#     18,
#     times,
#     folder,
#     "HETprofile",
#     "HET (μM N)",
#     colorrange = (0, 50),
#     colormap = Reverse(:afmhot),
# )


# record_vertical_tracer(
#     DOM,
#     z,
#     18,
#     times,
#     folder,
#     "DOMprofile",
#     "DOM (μM N)",
#     colorrange = (0, 50),
#     colormap = Reverse(:CMRmap),
# )

# record_vertical_tracer(
#     POM,
#     z,
#     18,
#     times,
#     folder,
#     "POMprofile",
#     "POM (μM N)",
#     colorrange = (0, 50),
#     colormap = Reverse(:greenbrownterrain),
# )

 record_vertical_tracer(
     NUT,
     z,
     18,
     times,
     folder,
     "NUTprofile",
     "Nutrients (μM N)",
     colorrange = (0, 25),
     colormap = Reverse(:cherry),
 )

 record_vertical_tracer(
     O₂,
     z,
     18,
     times,
     folder,
     "O2profile",
     "Dissolved Oxygen (μM)",
     colorrange = (0, 350),
     colormap = :turbo,
 )

record_vertical_tracer(
    O₂_relat,
    z,
    18,
    times,
    folder,
    "O₂_relat_profile",
    "O₂ relative (%)",
    colorrange = (0, 110),
    colormap = :gist_stern,
)

