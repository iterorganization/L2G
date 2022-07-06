"""This file holds the FieldLineTracer class which helps the user setup the
case.
"""
import l2g.comp.core
import l2g.comp
import l2g.utils.meshio
import l2g.equil
import l2g.hlm.general
import l2g.hlm.steady_state
import l2g.hlm.ramp_down
import numpy as np

import logging
from time import perf_counter
log = logging.getLogger(__name__)

from typing import Union, List

class FieldLineTracer:
    """FLT object for performing FLT on a given input mesh.

    Attributes:
        name (str): Name of the case.
        parameters (l2g.comp.Parameters): Holds parameters for the tracer.
        options (l2g.comp.Options): Options on what to run.
        flt_obj (l2g.comp.core.PyFLT): Cython wrapper over C++ kernel.
        embree_obj (l2g.comp.core.PyEmbreeAccell): Cython wrapper over C++
            Embree.
        equilibrium (l2g.equil.Equilibrium): Data class for holding equilibrium
            data.
        eq (l2g.equil.EQ): Equilibrium diagnostics class.
        point_results (l2g.comp.L2GPointsResults): Data class for FLT results
            performed on points (e.g., for owl connection length graph).
        mesh_results (l2g.comp.L2GResults): Data class for FLT results
            performed on meshs (e.g., main storage FLT results)
        fl_results (l2g.comp.L2GFLs): Data class for fieldline points.
        owl_conlen_data (np.ndarray): Outer wall connection length graph array
        target_vertices (np.ndarray | list): Location of points
        target_triangles (np.ndarray | list): Mesh triangles. List of IDs.
        target_points (np.ndarray | list): List of R/Z/Phi points.
        fl_ids (list): List of triangle elements from which we wish to obtain
            the trajectories of the field lines.
        rd_results (l2g.comp.L2GRampDownHLM): Data class for storing HLM
            results for Ramp-Down scenario
        ss_results (l2g.comp.L2GsteadyStateHLM): Data class for storing HLM
            results for Steady-State scenario
        sup_results (l2g.comp.L2GStartUpHLM): Data class for storing HLM
            results for Start-Up scenario

    """
    def __init__(self):
        self.name = "L2G_flt_case-1"

        self.hlm_params: l2g.comp.HLM                 = l2g.comp.HLM()
        self.parameters: l2g.comp.Parameters          = l2g.comp.Parameters()
        self.options:    l2g.comp.Options             = l2g.comp.Options()
        self.flt_obj:    l2g.comp.core.PyFLT          = l2g.comp.core.PyFLT()
        self.embree_obj: l2g.comp.core.PyEmbreeAccell = l2g.comp.core.PyEmbreeAccell()
        self.equilibrium:l2g.equil.Equilibrium        = l2g.equil.Equilibrium()
        self.eq:         l2g.equil.EQ                 = l2g.equil.EQ()

        self.eq._verbose = True # For additional output

        #  Results objects
        self.point_results: l2g.comp.L2GPointResults = l2g.comp.L2GPointResults() # Holds FLT result on points
        self.mesh_results:  l2g.comp.L2GResults =      l2g.comp.L2GResults() # Holds FLT result on mesh
        self.fl_results:    l2g.comp.L2GFLs =          l2g.comp.L2GFLs() # Holds fieldlines

        self.owl_conlen_data: np.ndarray = None

        # Geometry data for the target.
        self.target_vertices:  list = []
        self.target_triangles: list = []
        self.target_points:    list = []
        self.fl_ids:           list = []

        # HLM results
        self.hlm_results: l2g.comp.L2GResultsHLM = l2g.comp.L2GResultsHLM()

        # # RampDown results
        # self.rd_results: l2g.comp.L2GRampDownHLM = l2g.comp.L2GRamp
        # # Flat-top, Steady-State results
        # self.ss_results: l2g.comp.L2GSteadyStateHLM = l2g.comp.


    def setParameters(self, parameters: 'l2g.comp.Parameters') -> None:
        self.parameters = parameters

    def setFltObj(self, flt_obj: 'l2g.comp.core.PyFLT') -> None:
        self.flt_obj = flt_obj

    def setEmbreeObj(self, embree_obj: 'l2g.comp.core.PyEmbreeAccell') -> None:
        self.embree_obj = embree_obj
        if self.flt_obj:
            self.flt_obj.applyRT(self.embree_obj)

    def commitMeshesToEmbree(self, mesh_files: Union[List[str], str], dim_mul=1e-3) -> list:
        """Commits mesh from files to embree.

        Returns:
            ok (bool): Ok if finished

        Arguments:
            mesh_files: Either a single file or multiple files
            dim_mul (float): Dimension multiplier to convert mesh dimension to
                meters. Default 1e-3.
        """
        log.info("Commiting meshes to Embree object...")
        if isinstance(mesh_files, str):
            mesh_files = [mesh_files]

        out_ids = []

        for file in mesh_files:
            log.info(f"Commiting {file} to Embree object.")
            v, t = l2g.utils.meshio.readMesh(file)
            out_ids.append(self.embree_obj.commitMesh(v * dim_mul, t))
            # Dimension is in meters

        return out_ids

    def setTargetData(self, vertices:list = [], cells:list = []) -> None:
        """Sets new target data via function. In this case a target mesh is
        provided to the case via providing a list of nodes and a list of cells.

        With this the past mesh results are cleared.
        """
        self.target_vertices = vertices
        self.target_triangles = cells
        self.mesh_results.reset()

    def setTargetPoints(self, points: list = []) -> None:
        """Sets new target data via function. In this case target points are
        provided. This is used when we wish to analyze FLs on a ROI (i.e.,
        midplane).

        With this past point results are cleared.
        """
        self.target_points = points
        self.point_results.reset()

    def setEquilibrium(self, equilibrium: l2g.equil.Equilibrium) -> None:
        """Set the equilibrium data, which is propagated to the eq analyze
        class and the external FLT C++ code.
        """
        self.equilibrium = equilibrium
        self.eq.setEquilibrium(equilibrium)
        self.parameters.plasma_r_displ = 0.0
        self.parameters.plasma_z_displ = 0.0

    def applyParameters(self) -> None:
        """Propagates the parameters to the external FLT C++ code. Run this
        before running any FLT run functions.
        """
        log.info("Commiting parameters")
        # Set the Embree pointer
        self.flt_obj.applyRT(self.embree_obj)

        # Apply RKF45 accuracy settings
        self.flt_obj.setAbsError(self.parameters.abs_error)
        self.flt_obj.setRelError(self.parameters.rel_error)

        # Apply plasma shift
        plasma_r_displ, plasma_z_displ = self.parameters.plasma_r_displ, self.parameters.plasma_z_displ
        self.flt_obj.setShift(plasma_r_displ, plasma_z_displ)
        # Apply plasma shift to the EQ object
        self.eq.setDisplacement(plasma_r_displ, plasma_z_displ)
        # Apply toroidal angle settings.
        self.flt_obj.setTimeSpan(self.parameters.time_end,
                                self.parameters.time_step)
        # Apply maximum connection length settings
        self.flt_obj.setMaximumConnectionLength(
            self.parameters.max_connection_length)
        # Apply the self intersection avoidance length.
        self.flt_obj.setSelfIntersectionAvoidanceLength(
            self.parameters.self_intersection_avoidance_length)
        # Set FLT option
        self.flt_obj.setFltOption(self.options.switch_runFLT)

    def loadEq(self) -> None:
        """Loads the equilibrium data to the external FLT C++ code. Data is
        used for FL equations. Run this before running any FLT run functions.
        """
        log.info("Dumping equilibrium data to FLT object.")
        self.flt_obj.setNDIM(self.equilibrium.grid_dim_r,
                             self.equilibrium.grid_dim_z)
        self.flt_obj.setRARR(self.equilibrium.grid_r)
        self.flt_obj.setZARR(self.equilibrium.grid_z)
        self.flt_obj.setPSI(self.equilibrium.psi.flatten())
        self.flt_obj.setVacuumFPOL(self.equilibrium.fpol_vacuum)
        self.flt_obj.setFARR(self.equilibrium.fpol_flux)
        self.flt_obj.setFPOL(self.equilibrium.fpol)

        # self.flt_obj.setRARR(self.eq._R)
        # self.flt_obj.setZARR(self.eq._Z)
        # self.flt_obj.setPSI(self.eq._psi.flatten())
        # self.flt_obj.setFARR(self.eq._psi_fpol)
        # self.flt_obj.setFPOL(self.eq._fpol)
        # self.flt_obj.setVacuumFPOL(self.eq._fpol[-1])

        self.flt_obj.prepare() # Prepares the interpolators


    def evaluateEq(self) -> None:
        """Tries to automatically determine the type of the equilibrium
        """
        log.info("Evaluating equilibrium")
        self.eq.evaluate()
        log.info(f"Type: {self.eq.type_}, LCFS flux: {self.eq.psiLCFS} Webb/rad")

    def processDataOnMesh(self) -> None:
        """Function that populates the self.mesh_results variable with data
        required to run FLT.
        """
        if not len(self.target_vertices):
            log.error("No target vertices loaded! Stopping")
            return

        start = perf_counter()
        log.info("Processing data on mesh...")
        l2g.comp.core.processData(flt_obj=self.flt_obj,
            tg_vertices=self.target_vertices,
            tg_cells=self.target_triangles, results=self.mesh_results,
            dim_mul=self.parameters.target_dim_mul)
        log.info(f"Processing done in {perf_counter() - start:.2f} seconds")

    def processDataOnPoints(self) -> None:
        """Function that populates the self.points_results variable with data
        required to run FLT.
        """
        if not len(self.target_points):
            log.error("No target points loaded! Stopping")
            return

        start = perf_counter()
        log.info("Processing data on points...")
        l2g.comp.core.processDataOnPoints(flt_obj=self.flt_obj,
            tg_vertices=self.target_points, results=self.point_results,
            dim_mul=self.parameters.target_dim_mul)
        log.info(f"Processing done in {perf_counter() - start:.2f} seconds")

    def getFL(self) -> None:
        """Starts FL trace and creates result container.
        """

        if not len(self.fl_ids):
            log.error("No target points to trace FL from. Stopping")
            return

        if self.mesh_results.empty:
            # We need the BaryCenters.
            self.processDataOnMesh()

        # Just to be sure, convert the fl ids to numpy unsigned int 32
        self.fl_ids = np.asarray(self.fl_ids, np.uint32)

        # Make a check on the IDs, to see if the values are within interval
        # of triangles
        N_triangles = len(self.target_triangles)
        _stop = False
        for el in self.fl_ids:
            if el > N_triangles or el < 0:
                log.info(f"Skipping element id {el} as it points to no triangle")
                log.info("Please see the settings and remove invalid FL IDs.")
                _stop = True
        if _stop:
            return

        log.info("Getting FLs...")
        start = perf_counter()
        self.fl_results = l2g.comp.core.getFL(self.flt_obj, self.target_vertices,
            self.target_triangles, self.fl_ids, self.mesh_results,
            self.parameters.target_dim_mul, self.options.switch_getFL_with_FLT)
        log.info(f"Finished getting FLs in {perf_counter() - start} seconds.")

    def getFLOnPoint(self, R: float, Z: float, Theta: float) -> list:
        """Obtain FL points that goes through the input parameters R, Z, Theta.
        """
        points = l2g.comp.core.getFLOnPoint(self.flt_obj, R, Z, Theta,
                                            self.options.switch_getFL_with_FLT)

        # Now let's order them.
        points = points[0] + points[1][::-1]
        return points

    def runFltOnMesh(self) -> None:
        if not len(self.target_vertices):
            log.error("No target vertices loaded! Stopping")
            return

        if self.mesh_results.empty:
            self.processDataOnMesh()

        # Determine what kind of FLT. By toroidal angle or max conlen.

        log.info("Starting FLT on mesh...")
        start = perf_counter()
        l2g.comp.core.runFLT(self.flt_obj, self.target_vertices,
            self.target_triangles, self.parameters.num_of_threads,
            self.mesh_results, self.parameters.target_dim_mul)

        # Additionally create the mask which shows the final plasma pattern
        self.mesh_results.mask = self.applyShadowMask(
            np.ones(self.mesh_results.conlen.shape))

        log.info(f"Finished. It took {perf_counter() - start} seconds.")

    def runFltOnPoints(self) -> None:
        if not len(self.target_points):
            log.error("No target points loaded! Stopping")
            return
        if self.point_results.empty:
            self.processDataOnPoints()

        log.info("Starting FLT on points...")
        start = perf_counter()
        l2g.comp.core.runFLTonPoints(flt_obj=self.flt_obj,
            tg_vertices=self.target_points,
            user_num_threads=self.parameters.num_of_threads,
            results=self.point_results,
            dim_mul=self.parameters.target_dim_mul)

        log.info(f"Finished. It took {perf_counter() - start} seconds.")

    def obtainOwlConlenGraph(self) -> None:
        """Obtains the outer wall connection length graph. Basically following
        FLs on the midplane.
        """

        radial_points = 2000
        toroidal_points = 360

        points = self.createMidplanePoints(radialPoints=radial_points,
                                           toroidalSegments=toroidal_points)
        self.target_points = points
        # Reset result on points.
        self.point_results.reset()

        R_points = points[:3*2000:3]

        log.info("Obtaining OWL connection length graph")
        # Getting the outer wall midplane profile.
        self.runFltOnPoints()

        radial_points = 2000
        toroidal_points = 360

        psiLcfs = self.eq.psiLCFS
        Rb, _, _, Bpm = self.eq.getOWL_midplane()
        flux = self.point_results.flux[:2000]
        # Calculate drsep
        drsep = np.abs((flux - psiLcfs) / (Rb * Bpm))

        # Actual drsep. take from R_points
        drsep_r = R_points - R_points[0]

        conlenUp = self.point_results.conlenUp.reshape((toroidal_points,
                                                        radial_points))
        conlenDown = self.point_results.conlenDown.reshape((toroidal_points,
                                                            radial_points))
        # Columns are the points, laid in the poloidal direction.
        conlenUp_mean = np.mean(conlenUp, 0)
        conlenDown_mean = np.mean(conlenDown, 0)
        out = np.empty((radial_points, 4), np.float64)
        out[:, 0] = drsep
        out[:, 1] = drsep_r
        out[:, 2] = conlenUp_mean
        out[:, 3] = conlenDown_mean
        self.owl_conlen_data = out

    def alignGeometryWithLCFS(self) -> None:
        """In most cases the wall silhouette used in generating the equilibrium
        do not coincide with the geometry used in FLT. Therefore in cases, such
        as limiter cases, where every millimeter counts, we require to align
        the input geometry, so that the equilibrium LCFS is limiting also on
        the geometry.

       .. note::

          In order to use it, the diagnostic for the equilibrium must have a
          the points of the LCFS or boundary contour prepared so that there
          are points to compare.

        The LCFS contour is clipped to the closest element on the mesh.
        Afterwards the distance is set as the displacement for the radial and
        vertical direction. If there are any original displacements, let's say
        vertical already applied, this is already taken into account since the
        LCFS points are following the shift, therefore the underlying C++ code
        will receive the full correct displacement.

        .. todo::

           Find a good algorithm to obtain LCFS contour points with a good
           resolution. But this has to be done in EQ class.

        """
        log.info("Checking whether the input geometry is aligned with the" +
                 " contact point of the data")

        # Simplest way to check is to check the Psi values. Note that this
        # will not work great, if for instance the poloidal flux map contains
        # a series of contour "islands" where the values overlap, meaning we
        # have a set of magnetic surfaces that have the same poloidal flux map.
        if self.mesh_results.empty:
            log.error("No FLT data.")
            return

        # Calculate the initial drsep.
        log.info("Initial check for drsep.")
        self.calculateDrsep()

        # Get the element that is closest.
        el_min = np.argmin(self.mesh_results.drsep)
        dr_min = self.mesh_results.drsep[el_min]
        elr_min = self.mesh_results.baryCent[3 * el_min]
        elz_min = self.mesh_results.baryCent[3 * el_min + 1]
        # psi_min = self.mesh_results.flux[el_min]
        log.info(f"Closest element {el_min}: R={elr_min}")
        log.info(f"Closest distance from boundary: {dr_min}")

        # Using the stored LCFS contour in the diagnostic class we find the
        # shortest path and use it as radial and vertical displacement.

        point = np.array([elr_min, elz_min])

        # Calculate the distance**2 from the closest element to the closest
        # point on the LCFS contour.

        displ = self.eq.alignLcfsToPoint(point)


        prev_r_displ = self.parameters.plasma_r_displ
        prev_z_displ = self.parameters.plasma_z_displ
        self.parameters.plasma_r_displ = displ[0]
        self.parameters.plasma_z_displ = displ[1]
        log.info(f"Setting new plasma displacement to:")
        log.info(f"R={displ[0]}, Z={displ[1]}")
        # REAPPLY parameters
        self.applyParameters()
        # DO NOT CHANGE EQ data
        self.eq.setDisplacement(prev_r_displ, prev_z_displ)
        self.loadEq()
        self.processDataOnMesh()

    def calculateDrsep(self) -> None:
        """From the evaluated Flux data and the equilibrium data evaluate the
        radial distance along the midplane for the input target mesh data.
        """
        log.info("Calculating distance from the boundary on midplane for each mesh element.")

        # Evaluate the equilibrium
        self.eq.evaluate()
        # Get IWL or OWL parameters
        log.info(f"Calculating for side: {self.parameters.side}")
        if self.parameters.side == "iwl":
            Rb, _, _, Bpm = self.eq.getIWL_midplane()
        else: # owl
            Rb, _, _, Bpm = self.eq.getOWL_midplane()

        drsep = (self.mesh_results.flux - self.eq.psiLCFS) / (Rb * Bpm)

        # In case of some COCOS notations, the flux gradient direction can
        # go either away from the plasma or inside of the plasma. In other
        # words the the flux values are either multiplied with -1 or 1. Hence
        # the correction. Otherwise in case when the gradient goes inside the
        # drsep value will become negative as the flux values will be less then
        # the boundary value
        drsep *= self.equilibrium.psi_sign

        self.mesh_results.drsep = drsep

        if self.eq.type_ == "div":
            # Also evaluate the distance from the 2nd separatrix.
            if self.parameters.side == "iwl":
                Rb, _, _, Bpm = self.eq.getIWL_midplane(lcfs=self.eq.psiLCFS2)
            else: # owl
                Rb, _, _, Bpm = self.eq.getOWL_midplane(lcfs=self.eq.psiLCFS2)

            if self.eq.psiLCFS2:
                drsep2 = (self.mesh_results.flux - self.eq.psiLCFS2) / (Rb * Bpm)
                drsep2 *= self.equilibrium.psi_sign
                self.mesh_results.drsep2 = drsep2
            else:
                self.mesh_results.drsep2 = np.zeros(drsep.shape)


    def createMidplanePoints(self, which:str='owl', length:float=0.3,
            radialPoints:int=2000, toroidalSegments:int=360,
            angularStep:float=1) -> np.ndarray:
        """Based on the Equilibrium data, create points placed on the midplane
        plane, just outside the boundary contour (only first point should
        generally have connection lengths going to infinity, or free path on
        the boundary magnetic surface).

        Arguments:
            which (str): On which side, iwl = inner wall, owl = outer wall
                we wish to obtain the points
            length (float): Length in the major radius direction. That is,
                effectively to which drsep we wish to populate points. In
                meters.
            radialPoints (int): How many points in radial direction we wish to
                have.
            toroidalSegments (int): How many segments or lines of points do we
                wish to have in poloidal direction. The number of segments
                times the angularStep should at maximum be 360 degrees.
            angularStep (float): Angular distance between the segments. In
                degrees.

        Returns:
            points (array): 1D array of points: (R1, Z1, Phi1, R2, Z2, Phi2,...)
        """
        log.info("Creating points on the midplane")
        # First we need to obtain the R, Z point of the midplane, which sits
        # on the boundary magnetic surface.
        if self.eq is None:
            log.error("No EQ instance loaded! Stopping.")
            return

        self.eq.evaluate()

        # Obtain the R, Z point
        # Rb, Z, Btotal, Bpm = self.eq.get_midplane_info(which=which)
        Rb, Z, _, _ = self.eq.getOWL_midplane()
        log.info(f"Midplane height [m]: {Z}")
        log.info(f"Radial start [m]: {Rb}")
        log.info(f"Radial span [m]: {length}")
        log.info(f"Number of radial points per segment: {radialPoints}")
        log.info(f"Toroidal segments: {toroidalSegments}")
        log.info(f"Angular step: {angularStep}")
        log.info(f"Total angle: {angularStep * toroidalSegments}")

        if which == 'owl':
            # Points go in +R direction away from owl midplane
            radial_dir = 1
        else:
            # Points go in -R direction away from iwl midplane
            radial_dir = -1

        # Now let's start creating the points.

        total_num_of_points = radialPoints * toroidalSegments * 3

        # Note: for post-processing, re-shaping will come in handy, so think of
        # segments as rows and radial points as columns
        # points.reshape((toroidalSegments, radialPoints * 3))

        points = np.empty(total_num_of_points, np.float32)
        radial_step = length / (radialPoints - 1)
        phi = 0
        angularStep = np.deg2rad(angularStep) # Convert to radians
        # To ensure that the poloidal segments are rows
        for i in range(toroidalSegments):
            # Now to ensure radial points would be columns

            for j in range(radialPoints):
                r_point = Rb + radial_dir * j * radial_step

                # x = r_point * np.cos(phi)
                # y = r_point * np.sin(phi)
                offset = (i * radialPoints + j) * 3
                points[offset] = r_point
                points[offset + 1] = Z
                points[offset + 2] = phi

            phi += angularStep

        return points

    def applyHLM(self) -> None:
        """Applies an exponential plasma profile. Either single or double,
        depending on the type stored in self.hlm_params.
        """
        if self.equilibrium is None:
            log.error("No loaded equilibrium. Stopping")
            return

        if self.mesh_results.empty:
            log.error("No mesh results to apply Ramp Down on. Stopping")
            return

        if self.mesh_results.drsep is None:
            log.error("Mesh results are empty. Stopping.")
            return

        if not self.hlm_results.empty:
            self.hlm_results.reset()

        # Obtain the arrays from the FLT results
        drsep = self.mesh_results.drsep # In meters
        Bdot = np.abs(self.mesh_results.Bdot)
        n_cells = drsep.shape[0]
        bfield_mag = np.linalg.norm(self.mesh_results.BVec.reshape((n_cells, 3)), axis=1)
        self.eq.evaluate()

        # Get IWL or OWL parameters
        if self.parameters.side == "iwl":
            Rb, _, Btotal, Bpm = self.eq.getIWL_midplane()
        else: # owl
            Rb, _, Btotal, Bpm = self.eq.getOWL_midplane()

        if self.hlm_params.hlm_type == "single":
            q_par = l2g.hlm.general.single_exponential_psol(drsep=drsep,
                Bt=Btotal, Bpm=Bpm, Rb=Rb, P_sol=self.hlm_params.p_sol,
                F=0.5, lambda_q=self.hlm_params.lambda_q)
        elif self.hlm_params.hlm_type == "double":
            q_par = l2g.hlm.general.double_exponential_psol(drsep=drsep,
                Bt=Btotal, Bpm=Bpm, Rb=Rb,
                lambda_q_main=self.hlm_params.lambda_q_main,
                lambda_q_near=self.hlm_params.lambda_q_near,
                Rq=self.hlm_params.ratio, P_sol=self.hlm_params.p_sol, F=0.5)
        elif self.hlm_params.hlm_type == "custom":
            # We have points and profile
            q_par = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points, profile=self.hlm_params.profile)
        elif self.hlm_params.hlm_type == "elm":
            # The ELM data is loaded with an extra step.
            interELM = l2g.hlm.steady_state.inter_ELM(drsep, Rb, Btotal, Bpm,
                    Rb=self.hlm_params.r_break)
            elm = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points, profile=self.hlm_params.profile)
            q_par = interELM + elm

            self.hlm_results.additional_arrays.append(elm)
            self.hlm_results.additional_arrays.append(interELM)
        elif self.hlm_params.hlm_type == "ramp-down":
            # P_sol taken from IMAS, hopefully.
            Ip = self.equilibrium.Ip
            if Ip < self.hlm_params.ip_transition:
                lambda_q = l2g.hlm.ramp_down.decay_length_L_mode_diverted(
                    a = self.equilibrium.a, Ip=Ip,
                    Area=self.equilibrium.Area, R=self.equilibrium.mag_axis_r)
                q_par = l2g.hlm.general.single_exponential_psol(drsep, Btotal,
                    Bpm, Rb, lambda_q * 1e-3, self.equilibrium.Psol)
            else:
                lambda_q = float("NaN")
                q_par = np.zeros(drsep.shape)

            self.hlm_results.additional_arrays = [lambda_q]
        else:
            q_par = np.zeros(drsep.shape)
        expansion = bfield_mag / Btotal
        q_inc = q_par * Bdot / Btotal

        q_inc = self.applyShadowMask(q_par * Bdot / Btotal)

        self.hlm_results.q_inc = q_inc
        self.hlm_results.q_par = q_par
        self.hlm_results.flux_expansion = expansion

    def applyShadowMask(self, array: np.ndarray) -> np.ndarray:
        """This function applies the shadow mask to an input array.

        ONLY USED WITH MESH RESULTS!!!

        Several criteria are used for determining if a certain area is wetted:

         * Connection length of a field line:
           if conlen > cutoff = True else False
         * Geometries which mark fieldlines as shadowed. This is mainly used in
           regions where the FLs escape the chamber and are not captured by
           anything.

        Arguments:
            array (np.ndarray): A 1D array with N_cells elements on which the
                mask is applied

        Returns:
            masked_array (np.ndarray): A copy of the input array with the mask.

        """

        # First apply the connection length criteria

        out = None

        out = np.where(
            self.mesh_results.conlen >= self.parameters.cutoff_conlen,
            array, 0)

        # Now for every geom id in parameters.artificial_fl_catcher_geom_id
        # mark the fls as zero

        for geom_id in self.parameters.artificial_fl_catcher_geom_id:
            log.info(f"Masking all elements with geom_hit_ids == {geom_id} to zero.")
            out = np.where(self.mesh_results.geom_hit_ids != geom_id,
                           out, 0)
        return out