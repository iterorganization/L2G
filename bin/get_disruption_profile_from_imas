#!/usr/bin/python

description = """Script for obtaining the disruption profiles from the
Disruption IDS in the IMAS IMAS_DISRUPTIONS database.

The Disruption IDS has the profile of the power stored in

disruption.profiles_1d[:].power_density_conductive_losses

The only problem is that the units are W/m^3, therefore we only take the
silhouette of the function (normalize it) and scale it with the total power:

total_conductive_power = (disruption.global_quantities.power_ohm - \
    disruption.global_quantities.power_radiated_electrons_impurities)[0]

Additionally, time slices that have no defined profile (no power points on
magnetic surfaces outside the boundary), are ignored.
"""

import argparse

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-s', '--shot', metavar='SHOT', type=int,
                    help='Shot number',
                    required=True)
parser.add_argument('-r', '--run', metavar='RUN', type=int,
                    help='Run number', required=True)
parser.add_argument('-u', '--user', metavar='USER', type=str, default="public",
                    help='Username')
parser.add_argument('-d', '--device',
                    metavar="DEVICE", type=str, default="ITER_DISRUPTIONS", help='Device')
parser.add_argument('-o', '--output-directory', help="Output directory",
                    metavar="OUTPUT_DIRECTORY", default=".", )

args = parser.parse_args()

import os
import sys

SHOT=args.shot
RUN=args.run
USER=args.user
DEVICE=args.device
OUTPUT_DIRECTORY=os.path.abspath(os.path.expanduser(args.output_directory))

if not os.path.exists(OUTPUT_DIRECTORY):
    print(f"Output directory {OUTPUT_DIRECTORY} does not exist. Creating it")
    os.makedirs(OUTPUT_DIRECTORY)

import imas

entry = imas.DBEntry(shot=SHOT, run=RUN, user_name=USER, db_name=DEVICE,
                     backend_id=imas.imasdef.MDSPLUS_BACKEND)
code, _ = entry.open()

if code != 0:
    sys.exit(code)

summary = entry.get("summary")
wall = entry.get("wall")

DRSEPS = []
HEAT_LOADS = []
TIMES = []

import l2g.equil
import numpy as np

for i in range(summary.time.size):

    TIME = summary.time[i]
    disruption = entry.get_slice("disruption", TIME, imas.imasdef.CLOSEST_SAMPLE)
    equilibrium = entry.get_slice("equilibrium", TIME, imas.imasdef.CLOSEST_SAMPLE)

    total_conductive_power = (disruption.global_quantities.power_ohm - \
        disruption.global_quantities.power_radiated_electrons_impurities)[0]

    if total_conductive_power < 10.0:
        continue


    profile = disruption.profiles_1d[0].power_density_conductive_losses

    # Get the flux values of the profile
    disruption_grid_psi = disruption.profiles_1d[0].grid.psi[:] / (2*np.pi)
    disruption_profile = disruption.profiles_1d[0].power_density_conductive_losses

    # Convert the disruption_grid_psi to drsep
    # Create the profile, but first generate the drsep!
    eq_data = l2g.equil.getEquilibriumFromIMAS(equilibrium.time_slice[0],
                        equilibrium.vacuum_toroidal_field, wall, summary)

    eq = l2g.equil.EQ(eq_data)
    eq.evaluate()
    if eq.psiLCFS is None:
        print(f"Ignoring time {TIME} since equilibrium couldn't be evaluated.")
        continue

    Rb, Z, Btotal, Bpm = eq.get_midplane_info(which="owl")
    disruption_grid_drsep = (eq._psi_center_sign * disruption_grid_psi - eq.psiLCFS) / (Rb * Bpm)

    # Positive drsep filter
    positive_value_drsep_mask = disruption_grid_drsep >= 0.0
    filtered_drsep = disruption_grid_drsep[positive_value_drsep_mask]
    filtered_profile = disruption_profile[positive_value_drsep_mask]

    if filtered_drsep.size == 0:
        continue
    TIMES.append(TIME)

    # Now we trim the zero value power densities at the tail, leaving only one
    zero_power_density_index = np.where(filtered_profile == 0.0)[0]
    if zero_power_density_index.size:
        filtered_drsep = filtered_drsep[:zero_power_density_index[0]+1]
        filtered_profile = filtered_profile[:zero_power_density_index[0]+1]


    # Normalize the power density profile
    steps = np.diff(filtered_drsep)
    midpoints = (filtered_profile[:-1] + filtered_profile[1:]) * 0.5
    norm_value = np.sum(midpoints * steps)
    normalized_profile = filtered_profile / norm_value

    # Scaling to get head load amplitude.
    # Divide by area which is the band of the profile on the drsep.
    area = np.pi * ((Rb + filtered_drsep[-1])**2 - Rb**2)
    # 0.5 comes from the fact that heat loads goes up and down or from "both"
    # midplanes
    scale = 0.5 * total_conductive_power / area

    scaled_profile = normalized_profile * scale

    DRSEPS.append(filtered_drsep)
    HEAT_LOADS.append(scaled_profile)

    fname = os.path.join(OUTPUT_DIRECTORY, f"shot_{SHOT}_run_{RUN}_time_{TIME:.4f}.profile")
    npy_arr = np.zeros((2, filtered_drsep.size))
    npy_arr[0, :] = filtered_drsep
    npy_arr[1, :] = scaled_profile
    np.savetxt(fname, npy_arr)
    print(f"Profile at {TIME} saved to {fname}.")

# Create a graphics
import matplotlib
def get_rjet_cycler(n: int) -> "matplotlib.colors.LinearSegmentedColormap":
    return plt.cycler("color", matplotlib.cm.jet_r(np.linspace(0, 1, n)))

import matplotlib.pyplot as plt
f, ax = plt.subplots()
ax.set_prop_cycle(get_rjet_cycler(len(HEAT_LOADS)))

for i in range(len(HEAT_LOADS)):
    ax.plot(DRSEPS[i], HEAT_LOADS[i])

ax.set_ylabel(r"$q_{\parallel}$ $[\frac{W}{m^2}]$")
ax.set_xlabel(r"$\Delta_{sep}$ - radial distance along the midplane [m]")
ax.grid()
ax.set_yscale("log")

norm = matplotlib.colors.Normalize(vmin=TIMES[0], vmax=TIMES[-1])
plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="jet_r"), ax=ax,
             label="Time [s]")

fname = os.path.join(OUTPUT_DIRECTORY, f"shot_{SHOT}_run_{RUN}_disruption_profile.png")
f.savefig(fname)