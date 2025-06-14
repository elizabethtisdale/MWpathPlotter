import matplotlib.pyplot as plt
import openpyxl
import pandas as pd
from pyproj import CRS, Transformer
import math
from matplotlib.widgets import Button
import tkinter as tk
from tkinter import filedialog, messagebox





def get_utm_epsg_code(lon):
    """Calculates the appropriate UTM EPSG code for a given longitude and latitude."""

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

    #Define Coordinate Reference Systems (CRS) (thanks ChatGPT)
    crs_geographic = CRS("EPSG:4269")  # NAD83
    epsg_projected = get_utm_epsg_code(lon) #to get the code for the UTM zone we are in
    crs_projected = CRS(f"EPSG:{epsg_projected}")

    #Objects that will be used to transform the coordinates between systems
    to_utm = Transformer.from_crs(crs_geographic, crs_projected, always_xy=True)
    from_utm = Transformer.from_crs(crs_projected, crs_geographic, always_xy=True)

    #Project the original coordinate to UTM (meters)
    easting, northing = to_utm.transform(lon, lat)

    # Calculate new points and transform them back to geographic coordinates
    points = {}
    points['EAST'] = easting + distance_meters, northing
    points['WEST'] = easting - distance_meters, northing
    points['NORTH'] = easting, northing + distance_meters
    points['SOUTH'] = easting, northing - distance_meters

    # Return dictionary with lat/lon values
    return (points['WEST'][0], points['EAST'][0]), (points['SOUTH'][1], points['NORTH'][1]), to_utm

def get_list_of_paths(df, to_utm):
    '''
    To get the list of start and end point of each path
    :param df: dataframe with data read directly from the filtering excel
    :return: list of tuples of oring and end of each path, and additional info for each
    '''

    paths = []
    additional_info = []

    for i, row in df.iterrows():

        trans_lon, trans_lat = to_utm.transform(row['Lon (NAD83)'], row['Lat (NAD83)'])
        recei_lon, recei_lat = to_utm.transform(row['Lon (NAD83).1'], row['Lat (NAD83).1'])

        transmitter_coordinate = (trans_lat, trans_lon)
        receiver_coordinate = (recei_lat, recei_lon)

        paths.append((transmitter_coordinate, receiver_coordinate))
        additional_info.append((row['Callsign'], row['Path No.']))


    return paths, additional_info


def plot_microwave_paths(paths, x_range, y_range, legend_info, site_coords, radius):
    """
    Plots 2D lines for each microwave path within a defined axis range.

    Args:
        paths (list): A list of tuples, where each tuple contains two points:
                      the start and end points. Ex: [((x1_s, y1_s), (x1_e, y1_e)), ...].
        x_range (tuple): A tuple with the minimum and maximum limits for the X-axis. (min_x, max_x).
        y_range (tuple): A tuple with the minimum and maximum limits for the Y-axis. (min_y, max_y).
    """
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

    fig.canvas.mpl_connect('pick_event', on_pick)

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

    fig.canvas.mpl_connect('motion_notify_event', on_hover)

    # --- "Toggle All" Button Setup ---
    button_ax = fig.add_axes([0.82, 0.1, 0.1, 0.025])  # Define button position
    toggle_button = Button(button_ax, 'Toggle All', hovercolor='0.9')

    # Click handler for the button
    def on_toggle_clicked(event):
        if not all_path_lines: return
        # Determine target state from the first path's visibility
        target_visible = not all_path_lines[0].get_visible()

        # Set visibility for all lines and legend items
        for line in all_path_lines:
            line.set_visible(target_visible)
        for handle in path_legend_handles:
            handle.set_alpha(1.0 if target_visible else 0.2)

        fig.canvas.draw_idle()

    toggle_button.on_clicked(on_toggle_clicked)
    fig._button = toggle_button  # Keep a reference to prevent garbage collection

    # --- Final Plot Configuration ---
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_title("Microwave Paths Map (UTM)")
    ax.grid(True)
    ax.set_aspect('equal', adjustable='box')

    plt.show()

################MAIN#########################


def browse_for_file():
    """Opens a file dialog to select an .xlsm file."""
    filepath = filedialog.askopenfilename(
        title="Select Filtered MW Excel File",
        filetypes=(("Excel Macro-Enabled Workbook", "*.xlsm"), ("All files", "*.*"))
    )
    if filepath:
        entry_path.delete(0, tk.END)  # Clear the current entry
        entry_path.insert(0, filepath)  # Insert the new path

def start_generation():
    """Gets the path from the entry box and runs the main script logic."""
    filepath = entry_path.get()
    if not filepath:
        messagebox.showwarning("Warning", "Please select an Excel file first.")
        return
    mw_path_df = pd.read_excel(filepath, sheet_name='mw_Results', skiprows=1, usecols='F:J, O:Q')
    wb = openpyxl.load_workbook(filepath)
    ws = wb['mw_Input']

    site_lat = ws['B2'].value
    site_long = ws['B3'].value

    path_radius = ws['B4'].value  # in miles

    # Obtain plot axis ranges

    x_range, y_range, to_utm = calculate_axis_ranges(site_lat, site_long, path_radius + 1)

    site_coords = to_utm.transform(site_long, site_lat)
    paths, additional_info = get_list_of_paths(mw_path_df, to_utm)

    plot_microwave_paths(paths, x_range, y_range, additional_info, site_coords, path_radius)

# Create the main window
root = tk.Tk()
root.title("Microwave Path Plotter")
root.geometry("500x150") # Set initial size

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
button_generate = tk.Button(main_frame, text="Generate Interactive Drawing", command=start_generation, bg="#c8e6c9")
button_generate.grid(row=1, column=0, columnspan=3, pady=20)

# Start the Tkinter event loop
root.mainloop()

pass