import matplotlib

matplotlib.use('TkAgg')  #Use the Tkinter backend
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import matplotlib.pyplot as plt
import openpyxl
import pandas as pd
import pyproj
import math
from matplotlib.widgets import Button
import tkinter as tk
from tkinter import filedialog, messagebox
import os
import geopandas as gpd
from shapely.geometry import LineString


import sys
import os

# This block sets the necessary environment variables for GDAL and PROJ
# to work correctly when the application is packaged with PyInstaller. (Thanks Gemini)
if getattr(sys, 'frozen', False):
    # If the application is run as a bundle, the PyInstaller bootloader
    # extends the sys module by a flag frozen=True and sets the app path.

    # Add the unpacked C-libraries folder to the system's PATH
    os.environ['PATH'] = sys._MEIPASS + os.pathsep + os.environ['PATH']

    # Set the PROJ and GDAL data paths
    os.environ['PROJ_LIB'] = os.path.join(sys._MEIPASS, 'proj')
    os.environ['GDAL_DATA'] = os.path.join(sys._MEIPASS, 'gdal')


def get_utm_epsg_code(lon):
    """Calculates the appropriate UTM EPSG code for a given longitude."""
    utm_zone = math.floor((lon + 180) / 6) + 1
    # For the US, the northern hemisphere codes (326xx) are used.
    # Build the EPSG code string. "326" is the prefix for WGS 84 / Northern Hemisphere.
    epsg_code = f"326{utm_zone:02d}"
    return epsg_code


def calculate_axis_ranges(lat, lon, distance_miles):
    """
    Calculates 4 new coordinates (N, S, E, W) at a given distance
    from an initial point anywhere in the USA.
    """
    distance_meters = distance_miles * 1609.344
    # Define Coordinate Reference Systems (CRS)
    crs_geographic = pyproj.CRS("EPSG:4269")  # NAD83
    epsg_projected = get_utm_epsg_code(lon)  # to get the code for the UTM zone we are in
    crs_projected = pyproj.CRS(f"EPSG:{epsg_projected}")

    # Objects that will be used to transform the coordinates between systems
    to_utm = pyproj.Transformer.from_crs(crs_geographic, crs_projected, always_xy=True)

    # Project the original coordinate to UTM (meters)
    easting, northing = to_utm.transform(lon, lat)

    # Calculate new points
    points = {}
    points['EAST'] = easting + distance_meters
    points['WEST'] = easting - distance_meters
    points['NORTH'] = northing + distance_meters
    points['SOUTH'] = northing - distance_meters

    # Return UTM ranges and the transformer object
    x_range_utm = (points['WEST'], points['EAST'])
    y_range_utm = (points['SOUTH'], points['NORTH'])
    return x_range_utm, y_range_utm, to_utm, epsg_projected


def get_list_of_paths(df, to_utm):
    """
    Gets the list of start and end points for each path.
    :param df: DataFrame with data read directly from the filtering excel.
    :return: A tuple containing a list of path coordinate tuples and a list of additional info tuples.
    """
    paths = []
    additional_info = []

    for i, row in df.iterrows():
        # Transform geographic coordinates to UTM for plotting
        trans_easting, trans_northing = to_utm.transform(row['Lon (NAD83)'], row['Lat (NAD83)'])
        recei_easting, recei_northing = to_utm.transform(row['Lon (NAD83).1'], row['Lat (NAD83).1'])

        transmitter_coordinate = (trans_northing, trans_easting)
        receiver_coordinate = (recei_northing, recei_easting)

        paths.append((transmitter_coordinate, receiver_coordinate))
        additional_info.append((row['Callsign'], row['Path No.']))

    return paths, additional_info


def export_to_shapefile(full_output_path, paths, attributes, epsg_code):
    """
    Exports a list of microwave paths to a Shapefile using GeoPandas.
    This version does not require arcpy or ArcGIS Pro.
    """
    print(f"Starting GeoPandas export for {len(paths)} selected paths...")

    try:
        attrs_df = pd.DataFrame(attributes, columns=['Callsign', 'Path_Number'])
        # Create LineString geometries from the UTM coordinates
        geometries = [LineString([(path[0][1], path[0][0]), (path[1][1], path[1][0])]) for path in paths]
        gdf = gpd.GeoDataFrame(attrs_df, geometry=geometries, crs=f"EPSG:{epsg_code}")
        gdf['Path_Number'] = gdf['Path_Number'].astype(str)
        gdf.to_file(full_output_path, driver='ESRI Shapefile')

        print(f"Success! Shapefile saved to: {full_output_path}")
        messagebox.showinfo("Export Successful", f"{len(paths)} paths exported to:\n{full_output_path}")

    except Exception as e:
        print(f"An error occurred during export: {e}")
        messagebox.showerror("Export Error", f"Could not create the Shapefile.\nError: {e}")


# --- MODIFIED FUNCTION ---
def open_plot_window(parent, paths, x_range, y_range, legend_info, site_coords, radius, input_filepath, epsg_code):
    """
    Opens a NEW Tkinter window to display the Matplotlib plot.
    """
    plot_window = tk.Toplevel(parent)
    plot_window.title("Microwave Paths Map")
    plot_window.geometry("1500x1000")

    fig, ax = plt.subplots(figsize=(15, 10))
    fig.subplots_adjust(left=0.1, right=0.75, bottom=0.1, top=0.9)

    # Plot site and radius circle
    radius_meters = radius * 1609.344  # mi to m
    ax.plot(site_coords[0], site_coords[1], 'ro', markersize=10, label='Site Location')
    circle = plt.Circle((site_coords[0], site_coords[1]), radius_meters,
                        color='blue', fill=False, linestyle='--')
    ax.add_artist(circle)
    radius_handle = plt.Line2D([], [], color='blue', linestyle='--', label=f'{radius} mi Radius')

    # Plot each path individually
    all_path_lines = []
    for i, path in enumerate(paths):
        callsign, path_no = legend_info[i]
        start_point, end_point = path
        x_coords = [start_point[1], end_point[1]]
        y_coords = [start_point[0], end_point[0]]

        line, = ax.plot(x_coords, y_coords, marker='o', markersize=3, linestyle='-',
                        label=f'{callsign} No.{path_no}')

        # Store unique ID on the line object for the hover tooltip
        line.set_gid(f"{callsign} No.{path_no}")
        all_path_lines.append(line)

    # --- Interactive Legend Setup ---
    handles, labels = ax.get_legend_handles_labels()
    handles.insert(1, radius_handle)
    leg = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), title="Paths (Click to toggle)")

    # Map legend items to their corresponding plot lines
    legend_map = {}
    path_legend_handles = leg.legend_handles[2:]  # Skip Site and Radius handles
    for handle, line in zip(path_legend_handles, all_path_lines):
        handle.set_picker(5)  # Enable clicking on legend items
        legend_map[handle] = line

    # Click handler for toggling individual paths via the legend
    def on_pick(event):
        legend_handle = event.artist
        if legend_handle not in legend_map: return

        line = legend_map[legend_handle]
        line.set_visible(not line.get_visible())
        legend_handle.set_alpha(1.0 if line.get_visible() else 0.2)
        fig.canvas.draw_idle()

    # --- Hover Tooltip Setup ---
    annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="yellow", ec="black", lw=1))
    annot.set_visible(False)

    # Mouse move handler to show path details on hover
    def on_hover(event):
        if not event.inaxes == ax: return

        is_annot_visible = False
        for line in all_path_lines:
            if line.get_visible() and line.contains(event)[0]:
                annot.set_text(line.get_gid())
                annot.set_visible(True)
                annot.xy = (event.xdata, event.ydata)
                is_annot_visible = True
                break
        if not is_annot_visible and annot.get_visible():
            annot.set_visible(False)
        fig.canvas.draw_idle()

    # --- "Toggle All" Button Setup ---
    button_ax_toggle = fig.add_axes([0.80, 0.15, 0.15, 0.04])  # Define button position
    toggle_button = Button(button_ax_toggle, 'Toggle All', hovercolor='0.9')

    def on_toggle_clicked(event):
        if not all_path_lines: return
        target_visible = not all_path_lines[0].get_visible()
        for line in all_path_lines:
            line.set_visible(target_visible)
        for handle in path_legend_handles:
            handle.set_alpha(1.0 if target_visible else 0.2)
        fig.canvas.draw_idle()

    toggle_button.on_clicked(on_toggle_clicked)
    fig._toggle_button = toggle_button  # Keep a reference to prevent garbage collection

    # --- "Export" Button Setup ---
    button_ax_export = fig.add_axes([0.80, 0.10, 0.15, 0.04])
    export_button = Button(button_ax_export, 'Export Selected to Layer', hovercolor='0.9')

    # To export to an Esri Shapefile
    def on_export_clicked(event):
        visible_paths = []
        visible_attributes = []
        for i, line in enumerate(all_path_lines):
            if line.get_visible():
                visible_paths.append(paths[i])
                visible_attributes.append(legend_info[i])

        if not visible_paths:
            messagebox.showwarning("Warning", "No paths are selected/visible to export.")
            return

        # Save to the same path as the input excel
        input_dir = os.path.dirname(input_filepath)
        output_shp_path = os.path.join(input_dir, "Filtered_MW_paths.shp")
        export_to_shapefile(output_shp_path, visible_paths, visible_attributes, epsg_code)

    export_button.on_clicked(on_export_clicked)
    fig._export_button = export_button  # Keep a reference

    # --- Final Plot Configuration ---
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_title("Microwave Paths Map (UTM)")
    ax.grid(True)
    ax.set_aspect('equal', adjustable='box')

    # --- Tkinter Integration ---
    canvas = FigureCanvasTkAgg(fig, master=plot_window)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Add the Matplotlib toolbar (zoom, save, etc.)
    toolbar = NavigationToolbar2Tk(canvas, plot_window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Connect events to the new canvas
    canvas.mpl_connect('pick_event', on_pick)
    canvas.mpl_connect('motion_notify_event', on_hover)

    # DO NOT USE plt.show()
    # plt.show()


# --- Main Interface (MAIN) ---

def browse_for_file():
    """Opens a file dialog to select an Excel file."""
    filepath = filedialog.askopenfilename(
        title="Select Filtered MW Excel File",
        filetypes=(("Excel Files", "*.xlsx;*.xlsm"), ("All files", "*.*"))
    )
    if filepath:
        entry_path.delete(0, tk.END)  # Clear the current entry
        entry_path.insert(0, filepath)  # Insert the new path


# --- MODIFIED FUNCTION ---
def start_generation():
    """Gets the path from the entry box and runs the main script logic."""
    filepath = entry_path.get()
    if not filepath:
        messagebox.showwarning("Warning", "Please select an Excel file first.")
        return

    try:
        # Add a try/except block to catch any errors during loading/processing
        mw_path_df = pd.read_excel(filepath, sheet_name='mw_Results', skiprows=1, usecols='F:J, O:Q')
        wb = openpyxl.load_workbook(filepath, data_only=True)
        ws = wb['mw_Input']

        site_lat = ws['B2'].value
        site_long = ws['B3'].value
        path_radius = ws['B4'].value  # in miles

        # Obtain plot axis ranges
        x_range, y_range, to_utm, epsg_code = calculate_axis_ranges(site_lat, site_long, path_radius + 1)
        site_coords_utm = to_utm.transform(site_long, site_lat)
        paths, additional_info = get_list_of_paths(mw_path_df, to_utm)

        # Call the new function that opens the plot window
        open_plot_window(root, paths, x_range, y_range, additional_info, site_coords_utm, path_radius,
                         input_filepath=filepath, epsg_code=epsg_code)

    except Exception as e:
        # If something fails, show a clear error message.
        messagebox.showerror("Error", f"An error occurred while processing the file:\n\n{e}")


# --- Tkinter GUI Setup ---
root = tk.Tk()
root.title("Microwave Path Plotter")
root.geometry("500x150")  # Set initial size

# Create a frame for better organization
main_frame = tk.Frame(root, padx=10, pady=10)
main_frame.pack(fill=tk.BOTH, expand=True)

# --- Row 1: Label, Entry, and Browse Button ---
label_file = tk.Label(main_frame, text="Filtered MW Excel:")
label_file.grid(row=0, column=0, padx=5, pady=10, sticky="w")

entry_path = tk.Entry(main_frame, width=40)
entry_path.grid(row=0, column=1, padx=5, pady=10, sticky="ew")

button_browse = tk.Button(main_frame, text="Browse...", command=browse_for_file)
button_browse.grid(row=0, column=2, padx=5, pady=10)

# Configure the grid to make the entry box expand with the window
main_frame.grid_columnconfigure(1, weight=1)

# --- Row 2: Generate Button ---
button_generate = tk.Button(main_frame, text="Generate Plot", command=start_generation)
button_generate.grid(row=1, column=0, columnspan=3, pady=20)

# Start the Tkinter event loop
root.mainloop()