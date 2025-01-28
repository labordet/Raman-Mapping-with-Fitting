import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")  # For embedded matplotlib in Tk
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import filedialog, messagebox, StringVar, IntVar, BooleanVar, Radiobutton, Toplevel
from ttkthemes import ThemedTk
from tkinter import ttk

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import differential_evolution

# --------------------------------------------------------------------------------
# ALS BASELINE CORRECTION
# --------------------------------------------------------------------------------
def als_baseline_correction(spectrum, mask, lam=1e6, p=0.0001, niter=10):
    """
    Perform ALS baseline correction on the given 1D spectrum.
      - spectrum: 1D array of intensities
      - mask: Boolean array of the same length as spectrum, True where the baseline must pass
        through the original data
      - lam, p: regularization parameters
      - niter: number of iterations
    Returns the fitted baseline (1D array).
    """
    L = len(spectrum)
    # second-difference matrix
    D = diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for _ in range(niter):
        W = diags(w, 0, shape=(L, L))
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * spectrum)
        # Force baseline to match original data at mask points
        z[mask] = spectrum[mask]
        # Update weights
        w = p * (spectrum > z) + (1 - p) * (spectrum <= z)
    return z

# --------------------------------------------------------------------------------
# SINGLE LORENTZIAN + FIT
# --------------------------------------------------------------------------------
def single_lorentzian(x, h, pos, width):
    """
    Lorentzian peak shape:
    h / [1 + ((x - pos) / (0.5*width))^2]
    """
    return h / (1.0 + ((x - pos) / (0.5 * width)) ** 2)

def fit_lorentzian_differential_evolution(x, y, pos_bounds, width_bounds):
    """
    Fit a single Lorentzian peak to (x,y) data using differential evolution.
      - x, y: arrays of data
      - pos_bounds: (pos_min, pos_max)
      - width_bounds: (width_min, width_max)
    The bounds for height are set from 0.001*max(y) to 1.5*max(y).
    Returns best-fit parameters [h, pos, width].
    """
    def cost(params):
        h, pos, w = params
        y_fit = single_lorentzian(x, h, pos, w)
        return np.sum((y - y_fit) ** 2)

    max_y = np.max(y)
    h_min = 0.001 * max_y
    h_max = 1.5   * max_y

    bounds = [
        (h_min, h_max),
        (pos_bounds[0], pos_bounds[1]),
        (width_bounds[0], width_bounds[1])
    ]
    result = differential_evolution(
        cost,
        bounds,
        strategy='best1bin',
        popsize=15,
        mutation=(0.5, 1.0),
        recombination=0.7,
        maxiter=1000,
        polish=True
    )
    best_params = result.x  # [h, pos, w]
    return best_params

# --------------------------------------------------------------------------------
# MAIN APPLICATION
# --------------------------------------------------------------------------------
class SpectraAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Spectra Analysis - Advanced Version")
        # Make the main window a bit wider
        self.root.geometry("1200x800")
        self.root.resizable(True, True)

        # Main paned window: left side = controls, right side = plot
        self.main_paned = ttk.Panedwindow(self.root, orient="horizontal")
        self.main_paned.pack(fill="both", expand=True)

        # Left frame for all parameter widgets
        self.left_frame = ttk.Frame(self.main_paned)
        self.main_paned.add(self.left_frame, weight=0)

        # Right frame for the average spectrum plot
        self.right_frame = ttk.Frame(self.main_paned)
        self.main_paned.add(self.right_frame, weight=1)

        # Variables to store paths
        self.input_file = None
        self.output_dir = None

        # Data placeholders
        self.x = None
        self.spectra = None
        self.avg_spectrum = None

        # --------------------- THEME SELECTION ---------------------
        theme_frame = ttk.Frame(self.left_frame)
        theme_frame.pack(fill="x", padx=5, pady=5)

        ttk.Label(theme_frame, text="Select Theme:").pack(side="left", padx=5)
        self.theme_var = StringVar(value="arc")  # default theme
        self.themes = sorted(self.root.get_themes())  # available themes in ttkthemes
        self.theme_combo = ttk.Combobox(
            theme_frame, textvariable=self.theme_var,
            values=self.themes, width=15, state="readonly"
        )
        self.theme_combo.pack(side="left")
        self.theme_combo.bind("<<ComboboxSelected>>", self.change_theme)

        # --------------------- FILE SELECTION FRAME ---------------------
        self.file_frame = ttk.LabelFrame(self.left_frame, text="1) File and Folder Selection")
        self.file_frame.pack(fill="x", padx=10, pady=5)

        # Select .txt file
        ttk.Label(self.file_frame, text="Select .txt File:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.file_label = ttk.Label(self.file_frame, text="No file selected", width=50)
        self.file_label.grid(row=0, column=1, sticky="w")
        ttk.Button(self.file_frame, text="Browse", command=self.select_file).grid(row=0, column=2, padx=5)

        # Select output directory
        ttk.Label(self.file_frame, text="Select Output Directory:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.dir_label = ttk.Label(self.file_frame, text="No directory selected", width=50)
        self.dir_label.grid(row=1, column=1, sticky="w")
        ttk.Button(self.file_frame, text="Browse", command=self.select_directory).grid(row=1, column=2, padx=5)

        # --------------------- POINTS/LINES INPUTS ---------------------
        self.params_frame = ttk.LabelFrame(self.left_frame, text="2) Basic Parameters")
        self.params_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(self.params_frame, text="Points per Line:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.points_per_line_var = StringVar(value="110")
        self.points_per_line_entry = ttk.Entry(self.params_frame, width=10, textvariable=self.points_per_line_var)
        self.points_per_line_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.params_frame, text="Lines per Image:").grid(row=0, column=2, sticky="e", padx=5, pady=5)
        self.lines_per_image_var = StringVar(value="90")
        self.lines_per_image_entry = ttk.Entry(self.params_frame, width=10, textvariable=self.lines_per_image_var)
        self.lines_per_image_entry.grid(row=0, column=3, padx=5, pady=5)

        # --------------------- PEAK #1 FRAME ---------------------
        self.peak1_frame = ttk.LabelFrame(self.left_frame, text="3) Peak #1 Definition")
        self.peak1_frame.pack(fill="x", padx=10, pady=5)

        # Region
        ttk.Label(self.peak1_frame, text="Region Start:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.peak1_start_var = StringVar(value="240")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_start_var, width=10)\
            .grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.peak1_frame, text="Region End:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.peak1_end_var = StringVar(value="254")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_end_var, width=10)\
            .grid(row=0, column=3, padx=5, pady=5)

        # ALS lam, p
        ttk.Label(self.peak1_frame, text="ALS lam:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.peak1_lam_var = StringVar(value="1e6")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_lam_var, width=10)\
            .grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.peak1_frame, text="ALS p:").grid(row=1, column=2, padx=5, pady=5, sticky="e")
        self.peak1_p_var = StringVar(value="0.0001")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_p_var, width=10)\
            .grid(row=1, column=3, padx=5, pady=5)

        # Mask intervals
        ttk.Label(self.peak1_frame, text="Mask intervals (e.g. 170-175,327-338):")\
            .grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.peak1_mask_var = StringVar(value="")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_mask_var, width=30)\
            .grid(row=2, column=1, columnspan=3, padx=5, pady=5, sticky="w")

        # Lorentzian fit bounds
        ttk.Label(self.peak1_frame, text="Position bounds [min, max]:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.peak1_posmin_var = StringVar(value="245")
        self.peak1_posmax_var = StringVar(value="250")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_posmin_var, width=10)\
            .grid(row=3, column=1, padx=5, pady=5)
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_posmax_var, width=10)\
            .grid(row=3, column=2, padx=5, pady=5)

        ttk.Label(self.peak1_frame, text="Width bounds [min, max]:")\
            .grid(row=3, column=3, padx=5, pady=5, sticky="e")
        self.peak1_wmin_var = StringVar(value="1")
        self.peak1_wmax_var = StringVar(value="20")
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_wmin_var, width=10)\
            .grid(row=3, column=4, padx=5, pady=5)
        ttk.Entry(self.peak1_frame, textvariable=self.peak1_wmax_var, width=10)\
            .grid(row=3, column=5, padx=5, pady=5)

        # Button to preview baseline
        ttk.Button(
            self.peak1_frame,
            text="Preview Baseline Correction (Peak 1)",
            command=lambda: self.preview_baseline(peak=1)
        ).grid(row=4, column=0, columnspan=6, pady=5)

        # --------------------- PEAK #2 FRAME ---------------------
        self.peak2_frame = ttk.LabelFrame(self.left_frame, text="4) Peak #2 Definition")
        self.peak2_frame.pack(fill="x", padx=10, pady=5)

        # Region
        ttk.Label(self.peak2_frame, text="Region Start:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.peak2_start_var = StringVar(value="327")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_start_var, width=10)\
            .grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.peak2_frame, text="Region End:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.peak2_end_var = StringVar(value="338")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_end_var, width=10)\
            .grid(row=0, column=3, padx=5, pady=5)

        # ALS lam, p
        ttk.Label(self.peak2_frame, text="ALS lam:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.peak2_lam_var = StringVar(value="1e6")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_lam_var, width=10)\
            .grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.peak2_frame, text="ALS p:").grid(row=1, column=2, padx=5, pady=5, sticky="e")
        self.peak2_p_var = StringVar(value="0.0001")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_p_var, width=10)\
            .grid(row=1, column=3, padx=5, pady=5)

        # Mask intervals
        ttk.Label(self.peak2_frame, text="Mask intervals (e.g. 170-175,327-338):")\
            .grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.peak2_mask_var = StringVar(value="")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_mask_var, width=30)\
            .grid(row=2, column=1, columnspan=3, padx=5, pady=5, sticky="w")

        # Lorentzian fit bounds
        ttk.Label(self.peak2_frame, text="Position bounds [min, max]:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.peak2_posmin_var = StringVar(value="330")
        self.peak2_posmax_var = StringVar(value="335")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_posmin_var, width=10)\
            .grid(row=3, column=1, padx=5, pady=5)
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_posmax_var, width=10)\
            .grid(row=3, column=2, padx=5, pady=5)

        ttk.Label(self.peak2_frame, text="Width bounds [min, max]:")\
            .grid(row=3, column=3, padx=5, pady=5, sticky="e")
        self.peak2_wmin_var = StringVar(value="1")
        self.peak2_wmax_var = StringVar(value="20")
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_wmin_var, width=10)\
            .grid(row=3, column=4, padx=5, pady=5)
        ttk.Entry(self.peak2_frame, textvariable=self.peak2_wmax_var, width=10)\
            .grid(row=3, column=5, padx=5, pady=5)

        # Button to preview baseline
        ttk.Button(
            self.peak2_frame,
            text="Preview Baseline Correction (Peak 2)",
            command=lambda: self.preview_baseline(peak=2)
        ).grid(row=4, column=0, columnspan=6, pady=5)

        # --------------------- SUBTRACTION CHOICE ---------------------
        self.subtraction_frame = ttk.LabelFrame(self.left_frame, text="5) Subtraction Choice")
        self.subtraction_frame.pack(fill="x", padx=10, pady=5)

        self.sub_choice = IntVar(value=1)  # 1 => Peak1 - Peak2, 2 => Peak2 - Peak1
        Radiobutton(self.subtraction_frame, text="Peak1 - Peak2", variable=self.sub_choice, value=1)\
            .grid(row=0, column=0, padx=5, pady=5)
        Radiobutton(self.subtraction_frame, text="Peak2 - Peak1", variable=self.sub_choice, value=2)\
            .grid(row=0, column=1, padx=5, pady=5)

        # --------------------- NEW CHECKBOX: SAVE FITTED SPECTRA ---------------------
        self.save_fitted_spectra_var = BooleanVar(value=False)
        ttk.Checkbutton(self.left_frame, text="Save Fitted Spectra", variable=self.save_fitted_spectra_var)\
            .pack(padx=10, pady=5, anchor="w")

        # --------------------- ANALYZE BUTTON + PROGRESS BAR ---------------------
        self.analysis_frame = ttk.Frame(self.left_frame)
        self.analysis_frame.pack(fill="x", padx=10, pady=10)

        self.analyze_button = ttk.Button(self.analysis_frame, text="Analyze Full Data", command=self.analyze)
        self.analyze_button.pack(side="top", pady=5)

        self.progress = ttk.Progressbar(self.analysis_frame, orient="horizontal", length=250, mode="determinate")
        self.progress.pack(side="top", pady=5)

        # --------------------- PLOT FRAME for AVERAGE SPECTRUM ---------------------
        self.plot_frame = ttk.LabelFrame(self.right_frame, text="Average Spectrum")
        self.plot_frame.pack(fill="both", expand=True, padx=10, pady=10)

    # ----------------------------------------------------------------------------
    # THEME CHANGE
    # ----------------------------------------------------------------------------
    def change_theme(self, event=None):
        chosen = self.theme_var.get()
        self.root.set_theme(chosen)

    # ----------------------------------------------------------------------------
    # FILE/DIR SELECTION
    # ----------------------------------------------------------------------------
    def select_file(self):
        path = filedialog.askopenfilename(
            title="Select the .txt file",
            filetypes=[("Text files", "*.txt")]
        )
        if path:
            self.input_file = path
            self.file_label.config(text=os.path.basename(path))
            self.load_and_plot_average()

    def select_directory(self):
        path = filedialog.askdirectory(title="Select the directory to save results")
        if path:
            self.output_dir = path
            self.dir_label.config(text=path)

    # ----------------------------------------------------------------------------
    # DATA LOADING AND AVERAGE PLOT
    # ----------------------------------------------------------------------------
    def load_data(self, file_path):
        """
        Reads a space-delimited text file. 
        First column = X-axis; subsequent columns = spectra.
        Returns x_axis, df_of_spectra
        """
        # Use sep='\s+' to replace delim_whitespace=True (avoids FutureWarning).
        data = pd.read_csv(file_path, sep='\s+', header=None)
        x_axis = data.iloc[:, 0].values
        spectra = data.iloc[:, 1:]
        return x_axis, spectra

    def load_and_plot_average(self):
        if not self.input_file:
            return

        for w in self.plot_frame.winfo_children():
            w.destroy()

        self.x, self.spectra = self.load_data(self.input_file)
        self.avg_spectrum = self.spectra.mean(axis=1).values

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(self.x, self.avg_spectrum, label='Average Spectrum')
        ax.set_xlabel('X')
        ax.set_ylabel('Intensity')
        ax.set_title('Average Spectrum')
        ax.legend()
        ax.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")
        canvas.get_tk_widget().pack(fill="both", expand=True)

    # ----------------------------------------------------------------------------
    # PREVIEW BASELINE CORRECTION (AVERAGED SPECTRUM)
    # ----------------------------------------------------------------------------
    def preview_baseline(self, peak=1):
        if self.spectra is None or self.avg_spectrum is None:
            messagebox.showerror("No Data", "Please select a valid .txt file first.")
            return

        # Grab user inputs for the selected peak
        if peak == 1:
            start = float(self.peak1_start_var.get())
            end   = float(self.peak1_end_var.get())
            lam   = float(self.peak1_lam_var.get())
            p     = float(self.peak1_p_var.get())
            mask_str = self.peak1_mask_var.get()
        else:
            start = float(self.peak2_start_var.get())
            end   = float(self.peak2_end_var.get())
            lam   = float(self.peak2_lam_var.get())
            p     = float(self.peak2_p_var.get())
            mask_str = self.peak2_mask_var.get()

        x_full = self.x
        y_full = self.avg_spectrum

        region_mask = (x_full >= start) & (x_full <= end)
        x_region = x_full[region_mask]
        y_region = y_full[region_mask]

        # Build mask array from user intervals
        mask_array = np.zeros_like(x_region, dtype=bool)
        intervals = self.parse_mask_intervals(mask_str)
        for (mstart, mend) in intervals:
            inside = (x_region >= mstart) & (x_region <= mend)
            mask_array[inside] = True

        # Compute baseline
        baseline = als_baseline_correction(y_region, mask_array, lam=lam, p=p, niter=10)

        # Plot in a new window
        win = Toplevel(self.root)
        win.title(f"Preview Baseline - Peak {peak}")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x_region, y_region, label='Original Spectrum')
        ax.plot(x_region, baseline, label='Baseline', color='red')
        ax.set_title(f"Baseline Preview - Peak {peak}")
        ax.legend()

        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")

    def parse_mask_intervals(self, mask_str):
        """
        Parse a string like '170-175, 327-338' into a list of (start, end) float pairs.
        """
        intervals = []
        if not mask_str.strip():
            return intervals
        parts = mask_str.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                try:
                    s, e = part.split('-')
                    s, e = float(s), float(e)
                    if s > e:
                        s, e = e, s
                    intervals.append((s, e))
                except:
                    pass
        return intervals

    # ----------------------------------------------------------------------------
    # MAIN ANALYSIS
    # ----------------------------------------------------------------------------
    def analyze(self):
        if not self.input_file or not self.output_dir:
            messagebox.showerror("Error", "Please select both input file and output directory.")
            return

        # Validate numeric inputs
        try:
            ppline = int(self.points_per_line_var.get())
            lpi = int(self.lines_per_image_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid input for Points/Line or Lines/Image.")
            return

        # Grab all peak1 parameters
        try:
            p1_start = float(self.peak1_start_var.get())
            p1_end   = float(self.peak1_end_var.get())
            p1_lam   = float(self.peak1_lam_var.get())
            p1_p     = float(self.peak1_p_var.get())
            p1_mask_str = self.peak1_mask_var.get()

            p1_posmin = float(self.peak1_posmin_var.get())
            p1_posmax = float(self.peak1_posmax_var.get())
            p1_wmin   = float(self.peak1_wmin_var.get())
            p1_wmax   = float(self.peak1_wmax_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric input in Peak #1 definition.")
            return

        # Grab all peak2 parameters
        try:
            p2_start = float(self.peak2_start_var.get())
            p2_end   = float(self.peak2_end_var.get())
            p2_lam   = float(self.peak2_lam_var.get())
            p2_p     = float(self.peak2_p_var.get())
            p2_mask_str = self.peak2_mask_var.get()

            p2_posmin = float(self.peak2_posmin_var.get())
            p2_posmax = float(self.peak2_posmax_var.get())
            p2_wmin   = float(self.peak2_wmin_var.get())
            p2_wmax   = float(self.peak2_wmax_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric input in Peak #2 definition.")
            return

        # We'll do 2 * #spectra fits total
        num_spectra = self.spectra.shape[1]
        total_fits = 2 * num_spectra

        self.progress["value"] = 0
        self.progress["maximum"] = total_fits

        # Prepare arrays to store the fit results
        h1_vals = np.zeros(num_spectra)
        w1_vals = np.zeros(num_spectra)
        c1_vals = np.zeros(num_spectra)

        h2_vals = np.zeros(num_spectra)
        w2_vals = np.zeros(num_spectra)
        c2_vals = np.zeros(num_spectra)

        # Pre-parse the mask intervals
        p1_intervals = self.parse_mask_intervals(p1_mask_str)
        p2_intervals = self.parse_mask_intervals(p2_mask_str)

        # If user wants fitted spectra saved, create a subfolder
        save_fitted = self.save_fitted_spectra_var.get()
        fitted_folder = None
        if save_fitted:
            fitted_folder = os.path.join(self.output_dir, "Fitted_Spectra")
            os.makedirs(fitted_folder, exist_ok=True)

        # Iterate over each column (spectrum)
        for i in range(num_spectra):
            y_full = self.spectra.iloc[:, i].values

            # ---------------- PEAK #1 ----------------
            mask_region1 = (self.x >= p1_start) & (self.x <= p1_end)
            x1 = self.x[mask_region1]
            y1 = y_full[mask_region1]

            # Build the mask array for peak1
            mask_array1 = np.zeros_like(x1, dtype=bool)
            for (ms, me) in p1_intervals:
                inside = (x1 >= ms) & (x1 <= me)
                mask_array1[inside] = True

            # ALS baseline for peak1
            baseline1 = als_baseline_correction(y1, mask_array1, lam=p1_lam, p=p1_p, niter=10)
            y1_corrected = y1 - baseline1

            # Fit Lorentzian for peak1
            best_h1, best_c1, best_w1 = fit_lorentzian_differential_evolution(
                x1, y1_corrected,
                pos_bounds=(p1_posmin, p1_posmax),
                width_bounds=(p1_wmin, p1_wmax)
            )
            h1_vals[i] = best_h1
            w1_vals[i] = best_w1
            c1_vals[i] = best_c1

            self.progress["value"] += 1
            self.progress.update()

            # Optionally save a figure of the fitted result (Peak 1)
            if save_fitted:
                fig_p1, ax_p1 = plt.subplots(figsize=(6,4))
                ax_p1.plot(x1, y1, 'b-', label='Original Spectrum')
                ax_p1.plot(x1, baseline1, 'r--', label='Baseline')
                # Construct the fitted Lorentzian in original scale
                y1_fit = baseline1 + single_lorentzian(x1, best_h1, best_c1, best_w1)
                ax_p1.plot(x1, y1_fit, 'g-', label='Lorentzian Fit')
                ax_p1.set_title(f"Peak1 Fit - Spectrum {i+1}")
                ax_p1.legend()
                fig_p1.savefig(os.path.join(fitted_folder, f"peak1_spectrum_{i+1:04d}.png"), dpi=150)
                plt.close(fig_p1)

            # ---------------- PEAK #2 ----------------
            mask_region2 = (self.x >= p2_start) & (self.x <= p2_end)
            x2 = self.x[mask_region2]
            y2 = y_full[mask_region2]

            # Build the mask array for peak2
            mask_array2 = np.zeros_like(x2, dtype=bool)
            for (ms, me) in p2_intervals:
                inside = (x2 >= ms) & (x2 <= me)
                mask_array2[inside] = True

            baseline2 = als_baseline_correction(y2, mask_array2, lam=p2_lam, p=p2_p, niter=10)
            y2_corrected = y2 - baseline2

            best_h2, best_c2, best_w2 = fit_lorentzian_differential_evolution(
                x2, y2_corrected,
                pos_bounds=(p2_posmin, p2_posmax),
                width_bounds=(p2_wmin, p2_wmax)
            )
            h2_vals[i] = best_h2
            w2_vals[i] = best_w2
            c2_vals[i] = best_c2

            self.progress["value"] += 1
            self.progress.update()

            # Optionally save a figure of the fitted result (Peak 2)
            if save_fitted:
                fig_p2, ax_p2 = plt.subplots(figsize=(6,4))
                ax_p2.plot(x2, y2, 'b-', label='Original Spectrum')
                ax_p2.plot(x2, baseline2, 'r--', label='Baseline')
                y2_fit = baseline2 + single_lorentzian(x2, best_h2, best_c2, best_w2)
                ax_p2.plot(x2, y2_fit, 'g-', label='Lorentzian Fit')
                ax_p2.set_title(f"Peak2 Fit - Spectrum {i+1}")
                ax_p2.legend()
                fig_p2.savefig(os.path.join(fitted_folder, f"peak2_spectrum_{i+1:04d}.png"), dpi=150)
                plt.close(fig_p2)

        # ---------------- CREATE SUBFOLDERS AND SAVE RESULTS ----------------
        avg_c1 = np.mean(c1_vals)
        avg_c2 = np.mean(c2_vals)

        folder_peak1 = os.path.join(self.output_dir, f"Peak1_{avg_c1:.2f}")
        folder_peak2 = os.path.join(self.output_dir, f"Peak2_{avg_c2:.2f}")
        folder_sub = os.path.join(self.output_dir, "Substraction")

        os.makedirs(folder_peak1, exist_ok=True)
        os.makedirs(folder_peak2, exist_ok=True)
        os.makedirs(folder_sub,  exist_ok=True)

        # Save results for peak1
        self.save_peak_results(folder_peak1, h1_vals, w1_vals, c1_vals, ppline=ppline, lpi=lpi, peak_label="Peak1")
        # Save results for peak2
        self.save_peak_results(folder_peak2, h2_vals, w2_vals, c2_vals, ppline=ppline, lpi=lpi, peak_label="Peak2")

        # ---------------- SUBTRACTION ----------------
        # 1 => Peak1 - Peak2
        # 2 => Peak2 - Peak1
        if self.sub_choice.get() == 1:
            h_diff = h1_vals - h2_vals
            w_diff = w1_vals - w2_vals
            c_diff = c1_vals - c2_vals
        else:
            h_diff = h2_vals - h1_vals
            w_diff = w2_vals - w1_vals
            c_diff = c2_vals - c1_vals

        self.save_diff_results(folder_sub, h_diff, w_diff, c_diff, ppline=ppline, lpi=lpi)

        messagebox.showinfo("Success", "Analysis complete! Results saved.")

    # ----------------------------------------------------------------------------
    # HELPER: SAVE PEAK RESULTS
    # ----------------------------------------------------------------------------
    def save_peak_results(self, folder, h_vals, w_vals, c_vals, ppline, lpi, peak_label):
        """
        Saves heatmaps and an Excel of the (height, width, center) results for the given peak.
        """
        num_spectra = len(h_vals)
        try:
            h_map = h_vals.reshape((lpi, ppline))
            w_map = w_vals.reshape((lpi, ppline))
            c_map = c_vals.reshape((lpi, ppline))
        except ValueError:
            # If reshape fails, fallback to 1D arrays
            h_map = h_vals.reshape((1, num_spectra))
            w_map = w_vals.reshape((1, num_spectra))
            c_map = c_vals.reshape((1, num_spectra))

        # Height
        plt.figure()
        plt.imshow(h_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Height')
        plt.title(f'{peak_label} - Height')
        plt.savefig(os.path.join(folder, f'{peak_label}_Height.png'), dpi=150)
        plt.close()

        # Width
        plt.figure()
        plt.imshow(w_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Width')
        plt.title(f'{peak_label} - Width')
        plt.savefig(os.path.join(folder, f'{peak_label}_Width.png'), dpi=150)
        plt.close()

        # Center
        plt.figure()
        plt.imshow(c_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Center')
        plt.title(f'{peak_label} - Center')
        plt.savefig(os.path.join(folder, f'{peak_label}_Center.png'), dpi=150)
        plt.close()

        df = pd.DataFrame({
            "Height": h_vals,
            "Width": w_vals,
            "Center": c_vals
        })
        df.to_excel(os.path.join(folder, f'{peak_label}_results.xlsx'), index=False)

    # ----------------------------------------------------------------------------
    # HELPER: SAVE DIFFERENCE RESULTS
    # ----------------------------------------------------------------------------
    def save_diff_results(self, folder, h_diff, w_diff, c_diff, ppline, lpi):
        """
        Saves difference heatmaps and an Excel of the (height, width, center) differences.
        """
        num_spectra = len(h_diff)
        try:
            h_map = h_diff.reshape((lpi, ppline))
            w_map = w_diff.reshape((lpi, ppline))
            c_map = c_diff.reshape((lpi, ppline))
        except ValueError:
            h_map = h_diff.reshape((1, num_spectra))
            w_map = w_diff.reshape((1, num_spectra))
            c_map = c_diff.reshape((1, num_spectra))

        # Height diff
        plt.figure()
        plt.imshow(h_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Height Diff')
        plt.title('Height Difference')
        plt.savefig(os.path.join(folder, 'Height_Diff.png'), dpi=150)
        plt.close()

        # Width diff
        plt.figure()
        plt.imshow(w_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Width Diff')
        plt.title('Width Difference')
        plt.savefig(os.path.join(folder, 'Width_Diff.png'), dpi=150)
        plt.close()

        # Center diff
        plt.figure()
        plt.imshow(c_map, cmap='plasma', aspect='auto')
        plt.colorbar(label='Center Diff')
        plt.title('Center Difference')
        plt.savefig(os.path.join(folder, 'Center_Diff.png'), dpi=150)
        plt.close()

        df = pd.DataFrame({
            "HeightDiff": h_diff,
            "WidthDiff": w_diff,
            "CenterDiff": c_diff
        })
        df.to_excel(os.path.join(folder, 'Difference_results.xlsx'), index=False)

# --------------------------------------------------------------------------------
# START
# --------------------------------------------------------------------------------
if __name__ == "__main__":
    # Use ThemedTk so that we can apply ttk themes easily
    root = ThemedTk(theme="arc")  # default
    app = SpectraAnalysisApp(root)
    root.mainloop()