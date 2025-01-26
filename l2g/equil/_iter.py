import glob
import os

import logging
log = logging.getLogger(__name__)

import math

def truncate(number: float, digits: int) -> float:
    """Truncate number decimal places to the number of digits.
    """
    stepper = 10.0 ** digits
    return math.trunc(number * stepper) / stepper

def getBackUpIMASWallIds(shot=116000, run=4, db_name="ITER_MD",
        user_name="public"):
    """Gets the default Wall IDS machine description from the ITER ITER_MD
    machine description database. 116000/4
    """

    import imas
    import imas.imasdef

    db_entry = imas.DBEntry(shot=shot, run=run, db_name=db_name,
            user_name=user_name, backend_id=imas.imasdef.MDSPLUS_BACKEND)
    db_entry.open()

    return db_entry.get("wall")

class EquilibriumIterator(object):
    """Iterator object that has equilibriums and acts as a iterator over them.
    """

    def __init__(self) -> None:
        self.type = None

        #: List for storing :py:class:`l2g.equil.Equilibrium` objects.
        self._equilibriums: list['Equilibrium'] = []

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
        self.truncate_digits = 6

    def correctHelicity(self, val):
        self._correct_helicity = val

    def loadEqdskEquilibriums(self, l: list =[]) -> None:
        """Loads EQDSK G files as equilibriums.

        Arguments:
            l (list): List of EQDSK G files. Can contain * for globbing
        """
        from l2g.equil import getEquilibriumFromEQDSKG, EQDSKIO

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
            if not eqdsk.successfullRead:
                log.error(f"Failed reading {os.path.basename(file)}")
                raise Exception
            equilibrium = getEquilibriumFromEQDSKG(eqdsk,
                correct_helicity=self._correct_helicity)
            self._equilibriums.append(equilibrium)
            self._times.append(i)

    def loadIMASEquilibriums(self, d: dict = {}) -> None:
        """Loads equilibriums from IMAS. The dictionary arguments required::

         * user
         * device
         * version

        And for specifying the times we need either:

         * times

        or a defined time interval with:
         * time_start
         * time_end

        Optionally you can control the total number of steps with

         * time_samples

        which takes this many time steps from the set of times from previous
        two options (taking into account the ending and starting element).

        Example:

        .. code-block:: python

           d = {
             "user": "user_name", # Default public
             "device": "device_name", # Default iter
             "version": "data_version", # Default 3
             # For time steps the priority goes as following
             # First directly a list of time steps. Order is not required
             # time steps can appear multiple times
             "times" : [0.0, 0.1, 0.02, 0.02] # No default value
             # Or specifying the start and the end intervals of time to extract
             "time_start": 0 # Starting time, no default
             "time_end": 1 # End time, no default
             # And you can control the number of extracted steps, for both
             # cases by specifying time_samples
             "time_samples": 10 # Integer. No default value
           }

        Arguments:
            d (dict): Dictionary containing the IMAS parameters.
        """
        import numpy as np
        from l2g.equil import getEquilibriumFromIMAS

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
            if isinstance(times, np.ndarray):
                pass
            elif not isinstance(d["times"], list):
                times = list(times)

        time_samples = None
        if "time_samples" in d:
            time_samples = d["time_samples"]

        import imas
        interpolation = imas.imasdef.CLOSEST_INTERP

        # New API
        self._ids = imas.DBEntry(backend_id=imas.imasdef.MDSPLUS_BACKEND,
            db_name=device, shot=shot, run=run, user_name=user,
            data_version=version)
        self._ids.open()

        try:
            self._ids_wall = self._ids.get("wall")
        except:
            # Get backup wall IDS
            self._ids_wall = getBackUpIMASWallIds()


        log.info(f"Opening and reading data from SHOT={shot}, RUN={run}, device={device}, username={user}")

        self._ids_summary = self._ids.get("summary")

        # Get times. If times were not specifically defined, use the time_start
        # and time_end limits to extract the times.
        if times is None:

            times_indexes = np.where(np.logical_and(
                time_start <= self._ids_summary.time,
                self._ids_summary.time <= time_end))[0]
            times = self._ids_summary.time[times_indexes]

        # Extract only time_samples number of time steps. If it is defined.
        if time_samples and time_samples < len(times):
            idx = np.round(np.linspace(0, len(times) - 1, time_samples)).astype(int)
            times = times[idx]

        # Extract the times

        log.info(f"In total {len(times)} slices.")

        for t in times:
            log.info(f"Loading time slice {t}")
            self._ids_equilibrium = self._ids.get_slice("equilibrium", t,
                interpolation)
            self._ids_summary = self._ids.get_slice("summary", t,
                interpolation)
            time_slice = self._ids_equilibrium.time_slice[0]
            vacuum_toroidal_field = self._ids_equilibrium.vacuum_toroidal_field

            equilibrium = getEquilibriumFromIMAS(time_slice,
                vacuum_toroidal_field, self._ids_wall,
                self._ids_summary, correct_helicty=self._correct_helicity)

            # if interpolation == 1:
            #     # Closest interpolation
            #     t = time_slice.time

            self._equilibriums.append(equilibrium)

            # Get the actual time of the slice

            # In case in equilibrium the time_slice is written as -9e+40 or
            # 9e+40
            if abs(time_slice.time) > 1e10:
                truncate_time = truncate(t, self.truncate_digits)
            else:
                truncate_time = truncate(time_slice.time, self.truncate_digits)

            # Truncate the time
            self._times.append(truncate_time)

        log.info(f"Using {len(self._times)} slices as equilibrium input")

        return None

    def __len__(self) -> int:
        return len(self._equilibriums)

    def __getitem__(self, i):
        return i, self._times[i], self._equilibriums[i]

    def applyWallSilhouetteShift(self, r_shift: float, z_shift: float) -> None:
        """Applies a singular shift to all of the wall silhouettes stored
        inside.

        Arguments:
            r_shift (float): Radial shift. In meters
            z_shift (float): Vertical shift. In meters
        """
        for equilibrium in self._equilibriums:
            # Modify the wall silhouette points
            equilibrium.wall_contour_r = [_ + r_shift for _ in equilibrium.wall_contour_r]
            equilibrium.wall_contour_z = [_ + z_shift for _ in equilibrium.wall_contour_z]

    def applyPlasmaShift(self, r_shift: float | list[float], z_shift: float | list[float]) -> bool:
        """Applies a single or set of shifts (radial and/or vertical) to the
        plasma equilibrium. The number of shift should be the same as the
        number of equilibriums inside the iterator object.

        Arguments:
            r_shift (float): Radial shift. In meters
            z_shift (float): Vertical shift. In meters

        Returns:
            status(bool): True if it is appled okay, false otherwise.

        """

        log.info("Applying shift to input plasma equilibrium data")

        n_equilibriums = len(self._equilibriums)

        ok = True
        if (isinstance(r_shift, list) and len(r_shift) != n_equilibriums) or (not isinstance(r_shift, list) and n_equilibriums > 1 and r_shift != 0.0):
            ok = False
            log.error('You have not provided enough r shift values for all' +
                      ' instances of plasma! ')

        if (isinstance(z_shift, list) and len(z_shift) != n_equilibriums) or (not isinstance(z_shift, list) and n_equilibriums > 1 and z_shift != 0.0):
            ok = False
            log.error('You have not provided enough z shift values for all' +
                      ' instances of plasma!')
        if not ok:
            return False


        if isinstance(r_shift, list):
            for i,equilibrium in enumerate(self._equilibriums):
                log.info(f"Applying shift_r={r_shift[i]}m")
                equilibrium.mag_axis_r += r_shift[i]
                equilibrium.grid_r += r_shift[i]
            # It must be float
        elif r_shift != 0.0:
            log.info(f"Applying shift_r={r_shift} m")
            for equilibrium in self._equilibriums:
                # Apply shift to the equilibrium
                equilibrium.mag_axis_r += r_shift
                equilibrium.grid_r += r_shift

        if isinstance(z_shift, list):
            for i,equilibrium in enumerate(self._equilibriums):
                log.info(f"Applying shift_z={z_shift[i]}m")
                equilibrium.mag_axis_z += z_shift[i]
                equilibrium.grid_z += z_shift[i]
        elif z_shift != 0.0:
            log.info(f"Applying shift_z={z_shift} m")
            for equilibrium in self._equilibriums:
                # Apply shift to the equilibrium
                equilibrium.mag_axis_z += z_shift
                equilibrium.grid_z += z_shift
        return True
