# Nuclear Fission Plants Project

This repository collects all assignments developed for the Nuclear Fission Plants course. Each assignment is contained in its own folder and addresses a specific thermohydraulic design problem.

---

## Assignment 1

### 1.1 — Natural Circulation Loop Design (`project1.py`)
Determine the minimum elevation head (height) required to sustain natural circulation in the ISC (Intermediate Safety Circuit) piping loop. For each standard pipe diameter from the ASME table, the script computes the hydraulic pressure balance between driving buoyancy force and frictional/minor losses, producing a *Height vs. Diameter* optimization chart.

### 1.2 — Heat Exchanger Thermohydraulic Analysis (`project2.py`)
Size and verify the heat exchangers of the ISC and PSC (Primary Safety Circuit) loops. The script calculates the global heat-transfer coefficient *U*, the log-mean temperature difference, distributed and localized pressure drops along the circuit, and tracks the fluid temperature evolution. Results are exported to `result_ISC.csv` and `result_PSC.csv`.

---

## Two-phase flow laboratory
This laboratory exercise focuses on the analysis of two-phase (air-water) flow in a vertical pipe. The main script, `lab_fission.py`, performs the following tasks:
- It reads experimental data from `tab_dat_flowpat.csv`, which contains measurements of pressure, temperature, and flow rates.
- It calculates key experimental parameters, including mass flow rates for air and water, quality, mass flux, and superficial velocities.
- It determines the experimental void fraction and total pressure drop from the measured data.
- It compares these experimental results with the predictions of several widely-used correlations for void fraction (Homogeneous, Zivi, Chisholm, CISE, Drift-Flux) and pressure drop (including friction models like Friedel).
- The final comparison, containing both experimental and theoretical values, is exported to `risultati_analisi.csv`.

---

## Structure
- `assignment1/`: Scripts and data for Assignment 1.
  - `project1.py`: Natural circulation loop design (1.1).
  - `project2.py`: Heat exchanger thermohydraulic analysis (1.2).
  - `diameter_table.txt`: ASME standard pipe dimensions.
  - `result_ISC.csv`, `result_PSC.csv`: Output results.
- `lab/`: Scripts and data for the two-phase flow laboratory.
  - `lab_fission.py`: Main script for data analysis and comparison with correlations.
  - `tab_dat_flowpat.csv`: Raw experimental data.
  - `risultati_analisi.csv`: Output results comparing experimental data with correlations.
- `requirements.txt`: Lists the dependencies required to run the project.
- `README.md`: This file, providing an overview of the project.

## Setup
1. Create a virtual environment:
   ```bash
   python -m venv .venv
   ```
2. Activate the virtual environment:
   - On Windows:
     ```bash
     .venv\Scripts\activate
     ```
   - On macOS/Linux:
     ```bash
     source .venv/bin/activate
     ```
3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
Run the scripts in the `assignment1/` folder to perform the analyses. For example:
```bash
python assignment1/project1.py
```

### Git Workflow
To manage your project with Git, follow these steps:
- **Synchronize files**: Ensure your local repository is up-to-date by pulling the latest changes.
- **Modify files**: Make the necessary changes to your files.
- **Save changes**: Save the modified files (Control + S).
- **Stage changes**: Use the `+` icon next to the modified files under the "Changes" section to stage them.
- **Commit changes**: From the commit dropdown, enter a descriptive commit message indicating what was modified.
- **Push changes**: Save the commit and push it to the remote repository.

## License
This project is for educational purposes only.
