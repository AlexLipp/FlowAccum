"""
An example of how to use the D8Accumulator class to calculate drainage area and 
extract channel networks from a D8 flow direction raster. 
"""

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from d8_accumulator import D8Accumulator, write_geotiff, write_geojson

# Initialize the accumulator
accumulator = D8Accumulator("d8_eg.nc")
# Create an array of cell areas
trsfm = accumulator.ds.GetGeoTransform()
dx, dy = trsfm[1], trsfm[5] * -1
cell_area = np.ones(accumulator.arr.shape) * dx * dy
# Calculate drainage area in the base units of the D8 raster
print("Calculating drainage area")
drainage_area = accumulator.accumulate(weights=cell_area)

threshold = 1e7  # m^2
# Extracting channel segments from the drainage area array according to a threshold
print(f"Extracting channel segments with drainage area > {threshold} m^2")
channels = accumulator.get_channel_segments(drainage_area, threshold)

# Write the results to file
write_geotiff("drainage_area.tif", drainage_area, accumulator.ds)
print(f"Wrote drainage area to drainage_area.tif)")
write_geojson("channels.geojson", channels)
print(f"Wrote channels to channels.geojson")

# Pick a random point in the domain
startx, starty = 329000, 767500
print(f"Finding the channel profile starting at ({startx},{starty})")
# Find the corresponding node ID
start_node = accumulator.coord_to_node(startx, starty)
# Get the profile of the channel starting at this node
profile, distance = accumulator.get_profile(start_node)
area_on_profile = drainage_area.flatten()[profile]
profile_coords = [accumulator.node_to_coord(node) for node in profile]

# Find the sink node of the channel network
sink_node = accumulator.get_sink(start_node)
sink_x, sink_y = accumulator.node_to_coord(sink_node)

# Find the drainage mask of a different catchment
startx_2, starty_2 = 355000, 737500
print(f"Finding the basin mask that contains the point ({startx_2},{starty_2})")
sink_node_2 = accumulator.get_sink(accumulator.coord_to_node(startx_2, starty_2))
sink_x_2, sink_y_2 = accumulator.node_to_coord(sink_node_2)
upstream = accumulator.get_upstream_nodes(sink_node_2)
mask = accumulator.get_node_mask(upstream)
# Convert the Falses into np.NaN for transparent visualization
mask = mask.astype(float)
mask[mask == 0] = np.NaN

print("Visualizing results")
plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
plt.imshow(drainage_area, norm=LogNorm(), extent=accumulator.extent, cmap="Greys")
cb = plt.colorbar()
cb.set_label("Drainage area ($m^2$)")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.title("Drainage area ($m^2$)")

plt.subplot(2, 2, 2)
plt.imshow(drainage_area, norm=LogNorm(), extent=accumulator.extent, cmap="Greys")
cb = plt.colorbar()
cb.set_label("Drainage area ($m^2$)")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
for channel in channels.coordinates:
    x, y = zip(*channel)
    plt.plot(x, y, c="red")
plt.title(f"Channels with drainage area > {threshold} $m^2$")

plt.subplot(2, 2, 3)
plt.imshow(drainage_area, norm=LogNorm(), extent=accumulator.extent, cmap="Greys")
cb = plt.colorbar()
cb.set_label("Drainage area ($m^2$)")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.plot(
    [x for x, _ in profile_coords],
    [y for _, y in profile_coords],
    "g-",
    label="Channel profile",
)
plt.plot([startx], [starty], "bo", label="Start point (profile)")
plt.plot([sink_x], [sink_y], "ro", label="Sink point (profile)")

# Plot the mask of the second basin
plt.imshow(mask, extent=accumulator.extent)
# Plot the start point and the sink point as squares
plt.plot(startx_2, starty_2, "bs", label="Start point (mask)")
plt.plot(sink_x_2, sink_y_2, "rs", label="Sink point (mask)")
plt.title("Channel planform + basin mask")
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(distance, area_on_profile, "k-")
plt.xlabel("Distance from mouth (m)")
plt.ylabel("Drainage area ($m^2$)")
# Set the y and x limits to be the min max of the profile
plt.xlim(min(distance), max(distance))
plt.ylim(min(area_on_profile), max(area_on_profile))
plt.title("Channel profile (vertical)")

plt.tight_layout()
plt.show()
