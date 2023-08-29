from l2g.equil import (getEquilibriumFromIMAS, getEquilibriumFromEQDSKG,
                       EQDSKIO, Equilibrium)
import glob
import os

import logging
log = logging.getLogger(__name__)

import math

from typing import List

def truncate(number: float, digits: int) -> float:
    """Truncate number decimal places to the number of digits.
    """
    stepper = 10.0 ** digits
    return math.trunc(number * stepper) / stepper

class EquilibriumIterator(object):
    """Iterator object that has equilibriums and acts as a iterator over them.
    """

    def __init__(self) -> None:
        self.type = None

        #: List for storing :pyclass:`l2g.equil.Equilibrium` objects.
        self._equilibriums: List[Equilibrium] = []

        #: List of associated times.
        self._times = []

        #: IMAS IDS reference
        self._ids = None

        #: Wall IDS reference
        self._wall_ids = None

        #: Equilibrium IDS reference
        self._wall_equilibrium = None

        #: Flag for correcting helicity of equilibriums.
        self._correct_helicity = True

        #: Number of digits of time to write into identifiers
        self.truncate_digits = 3

    def correctHelicity(self, val):
        self._correct_helicity = val

    def loadEqdskEquilibriums(self, l: list =[]) -> None:
        """Loads EQDSK G files as equilibriums.

        Arguments:
            l (list): List of EQDSK G files. Can contain * for globbing
        """
        if isinstance(l, str):
            l = [l]

        eqdsk_files = []
        for file in l:
            if "*" in file:
                eqdsk_files += glob.glob(file)
            else:
                eqdsk_files.append(file)

        for i, file in enumerate(eqdsk_files):
            log.info(f"Loading {os.path.basename(file)}")
            eqdsk = EQDSKIO(file)
            equilibrium = getEquilibriumFromEQDSKG(eqdsk,
                correct_helicity=self._correct_helicity)
            self._equilibriums.append(equilibrium)
            self._times.append(i)

    def loadIMASEquilibriums(self, d: dict = {}) -> None:
        """Loads equilibriums from IMAS.
        """
        shot = d['shot']
        run = d['run']

        if 'user' not in d:
            user = 'public'
        else:
            user = d['user']

        if 'device' not in d:
            device = 'iter'
        else:
            device = d['device']

        if 'version' not in d:
            version = '3'
        else:
            version = d['version']

        # if 'times' in d:
        #     time_slices = d['times']
        # else:
        #     n_steps = int((d['time_end'] - d['time_start']) / d['time_step']) + 1
        #     time_slices = np.linspace(d['time_start'], d['time_end'], n_steps)

        # Ignore times, time_step and focus on time_start and time_end
        time_start = None
        if "time_start" in d:
            time_start = d["time_start"]

        time_end = None
        if "time_end" in d:
            time_end = d["time_end"]

        times = None
        if "times" in d:
            times = d["times"]
            if not isinstance(d["times"], list):
                times = [times]

        # OLD API
        # self._ids = imas.ids(shot, run)
        # self._ids.open_env(user, device, version)

        # self._ids_wall = self._ids.wall
        # self._ids_wall.get()

        # self._ids_summary = self._ids.summary
        # self._ids_summary.get()

        # self._ids_equilibrium = self._ids.equilibrium

        # Trouble in some cases
        # interpolation = imas.imasdef.INTERPOLATION
        # Closest interpolation.
        # interpolation = imas.imasdef.CLOSEST_SAMPLE
        import imas
        interpolation = imas.imasdef.CLOSEST_INTERP
        import numpy as np

        # New API
        self._ids = imas.DBEntry(backend_id=imas.imasdef.MDSPLUS_BACKEND,
            db_name=device, shot=shot, run=run, user_name=user,
            data_version=version)
        self._ids.open()

        self._ids_wall = self._ids.get("wall")


        log.info(f"Opening and reading data from SHOT={shot}, RUN={run}, device={device}, username={user}")

        self._ids_summary = self._ids.get("summary")

        # Get times
        if times is None:

            times_indexes = np.where(np.logical_and(
                time_start <= self._ids_summary.time,
                self._ids_summary.time <= time_end))[0]
            times = self._ids_summary.time[times_indexes]

        # Extract the times

        log.info(f"In total {len(times)} slices.")

        for t in times:
            log.info(f"Loading time slice {t}")
            self._ids_equilibrium = self._ids.get_slice("equilibrium", t,
                interpolation)
            self._ids_summary = self._ids.get_slice("summary", t,
                interpolation)
            slice = self._ids_equilibrium.time_slice[0]

            equilibrium = getEquilibriumFromIMAS(slice, self._ids_wall,
                self._ids_summary, correct_helicty=self._correct_helicity)

            if interpolation == 1:
                # Closest interpolation
                t = slice.time

            self._equilibriums.append(equilibrium)

            truncate_time = truncate(t, self.truncate_digits)
            # Truncate the time
            self._times.append(truncate_time)

        log.info(f"Using {len(self._times)} slices as equilibrium input")

        return None

    def __len__(self) -> int:
        return len(self._equilibriums)

    def __getitem__(self, i):
        return i, self._times[i], self._equilibriums[i]

    def applyWallSilhouetteShift(self, r_shift: float, z_shift: float):
        for equilibrium in self._equilibriums:
            # Modify the wall silhouette points
            equilibrium.wall_contour_r = [_ + r_shift for _ in equilibrium.wall_contour_r]
            equilibrium.wall_contour_z = [_ + z_shift for _ in equilibrium.wall_contour_z]

    def applyPlasmaShift(self, r_shift: float | List[float], z_shift: float | List[float]):

        log.info("Applying shift to input plasma equilibrium data")
        if isinstance(r_shift, list):
            if not len(r_shift) == len(self):
                log.error('You have not provided enough shift values for all' +
                          ' instances of plasma!')
                return

            for i,equilibrium in enumerate(self._equilibriums):
                log.info(f"Applying shift_r={r_shift[i]}m shift_z={z_shift[i]}m")
                equilibrium.mag_axis_r += r_shift[i]
                equilibrium.mag_axis_z += z_shift[i]
                equilibrium.grid_r += r_shift[i]
                equilibrium.grid_z += z_shift[i]
        else:
            log.info(f"Applying shift_r={r_shift}m shift_z={z_shift}m")
            for equilibrium in self._equilibriums:
                # Apply shift to the equilibrium
                equilibrium.mag_axis_r += r_shift
                equilibrium.mag_axis_z += z_shift
                equilibrium.grid_r += r_shift
                equilibrium.grid_z += z_shift
