from l2g.equil import Equilibrium
import glob
import os
from typing import Any

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


    An example of loading equilibrium data to the object and using it:


    .. code-block::

       import l2g.equil

       eqit = l2g.equil.EquilibriumIterator()

       eqit.loadEqdskEquilibriums(list_of_eqdsk_file_paths)
       # or
       eqit.loadIMASEquilibriums(dict_of_imas_params)

       for index, time, equilibrium in eqit:
           # Do work here
           pass

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
        """Creates Equilibrium objects from  EQDSK G files and stores them.

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
            if not os.path.exists(file):
                log.error(f"File {file} does not exist.")
                raise Exception

            log.info(f"Loading {os.path.basename(file)}")
            eqdsk = EQDSKIO(file)
            if not eqdsk.successfullRead:
                log.error(f"Failed reading {os.path.basename(file)}")
                raise Exception
            equilibrium = getEquilibriumFromEQDSKG(eqdsk,
                correct_helicity=self._correct_helicity)
            self._equilibriums.append(equilibrium)
            self._times.append(i)

    def loadIMASEquilibriums(self, d: dict[str, Any]) -> None:
        """Generate Equilibrium objects from IMAS.

        The dictionary arguments
        required::

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
             "shot": 135011, # Must be specified. Int
             "run": 7, # Must be specified. Int
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

    def __getitem__(self, i: int) -> tuple[int, float, "Equilibrium"]:
        return i, self._times[i], self._equilibriums[i]

    def applyBtMultiplier(self, bt_multiplier: float | list[float]) -> None:
        """Applies a multiplier to the Bt (fpol_vacuum) value. If the
        multiplier is a list of values, the number of values must correspond
        to the number of equilibriums.

        Arguments:
            bt_multiplier (float | list[float]): Float or a list of floats.
        """
        n_equilibriums = len(self._equilibriums)
        is_list = isinstance(bt_multiplier, list)

        if (is_list and len(bt_multiplier) != n_equilibriums):
            log.error('You have not provided enough bt_multiplier values for all' +
                      ' instances of plasma! ')
            raise Exception("Not enough Bt multiplier values provided")

        if is_list:
            for i, equilibrium in enumerate(self._equilibriums):
                equilibrium.fpol_vacuum *= bt_multiplier[i]
        else:
            for equilibrium in self._equilibriums:
                equilibrium.fpol_vacuum *= bt_multiplier

    def applyPsiMultiplier(self, psi_multiplier: float | list[float]) -> None:
        """Applies a multiplier to the Psi (fpol_vacuum) value. If the
        multiplier is a list of values, the number of values must correspond
        to the number of equilibriums.

        Arguments:
            psi_multiplier (float): Float.
        """
        n_equilibriums = len(self._equilibriums)
        is_list = isinstance(psi_multiplier, list)

        if (is_list and len(psi_multiplier) != n_equilibriums):
            log.error('You have not provided enough psi_multiplier values for all' +
                      ' instances of plasma! ')
            raise Exception("Not enough Psi multiplier values provided")

        if is_list:
            for i, equilibrium in enumerate(self._equilibriums):
                equilibrium.psi *= psi_multiplier[i]
        else:
            for equilibrium in self._equilibriums:
                equilibrium.psi *= psi_multiplier

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

    def applyPlasmaShift(self, r_shift: float | list[float],
                         z_shift: float | list[float]) -> None:
        """Applies a single or set of shifts (radial and/or vertical) to all
        loaded plasma equilibrium. If the shifts are a list of values, the
        number of values must correspond to the number of equilibriums.

        Arguments:
            r_shift (float): Radial shift. In meters
            z_shift (float): Vertical shift. In meters

        Returns:
            status(bool): True if it is applied okay, false otherwise.

        """

        log.info("Applying shift to input plasma equilibrium data")

        n_equilibriums = len(self._equilibriums)

        is_r_list = isinstance(r_shift, list)
        is_z_list = isinstance(z_shift, list)
        if is_r_list and len(r_shift) != n_equilibriums:
            raise Exception("Not enough radial shift values provided.")

        if is_z_list and len(z_shift) != n_equilibriums:
            raise Exception("Not enough vertical shift values provided.")


        if is_r_list:
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

        if is_z_list:
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
        return
