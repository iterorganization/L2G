import sys
import os

import numpy as np

from typing import Optional

import logging
log = logging.getLogger(__name__)

class DATA_BLOCK:
    required_keys: list = []
    optional_keys: list = []

    def __init__(self):
        self.data = {}

    def check_data(self, data: dict, required_keys=None) -> None:
        if required_keys is None:
            required_keys = self.required_keys
        missing_keys = []
        for key in required_keys:
            if key not in data:
                missing_keys.append(key)

        if missing_keys:
            log.error("Missing data in geometry block: " +
                          f"{','.join(missing_keys)}")
            sys.exit(1)

    def set_required_attr(self, data: dict, required_keys=None) -> None:
        if required_keys is None:
            required_keys = self.required_keys

        for key in required_keys:
            self.data[key] = data[key]

    def set_optional_attr(self, data: dict, optional_keys=None) -> None:
        if optional_keys is None:
            optional_keys = self.optional_keys

        for key in optional_keys:
            if key in data:
                self.data[key] = data[key]

    def load_data(self, data: dict) -> None:
        self.check_data(data)
        self.set_required_attr(data)
        self.set_optional_attr(data)


class GEOMETRY(DATA_BLOCK):
    """Class for storing information on FLT to be used
    """
    required_keys = ["name", "target_mesh", "shadow_meshes"]
    optional_keys = ["cutoff", "parameters", "fl_ids", "exclude_meshes",
                     "include_target_in_shadow", "afl_catcher_meshes",
                     "rot_axes", "rot_theta", "align_lcfs"]

    def __init__(self):
        super(GEOMETRY, self).__init__()

        self.data: dict = {"fl_ids": []}

class EQUILIBRIUM(DATA_BLOCK):
    required_keys: list = ["name", "equilibrium_type"]
    optional_keys: list = ["eqdsk_files", "imas", "custom_wall_limiter",
                           "wall_limiter_r", "wall_limiter_z",
                           "wall_silh_r_shift", "wall_silh_z_shift",
                           "plasma_r_displ", "plasma_z_displ"]
    def __init__(self):
        super(EQUILIBRIUM, self).__init__()
        self.data: dict = {"custom_wall_limiter": False}

    def check_data(self, data: dict) ->None:
        super(EQUILIBRIUM, self).check_data(data)
        # Now check combinations.
        if data["equilibrium_type"] == "eqdsk_files":
            if "eqdsk_files" not in data:
                log.error("Missing entry eqdsk_files in equilibrium block.")
                sys.exit(2)
        else:
            if "imas" not in data:
                log.error("Missing imas in equilibrium block")
                sys.exit(2)

    def load_data(self, data: dict) -> None:
        self.check_data(data)
        self.set_required_attr(data)
        self.set_optional_attr(data)


class HLM(DATA_BLOCK):
    required_keys: list = ["name", "hlm_type"]
    optional_keys: list = ["ip_transition", "r_break", "parameters", "ratio"]
    elm_required: list = ["shadow_meshes"]
    single_required: list = ["p_sol", "lambda_q"]
    double_required: list = ["p_sol", "lambda_q_near", "lambda_q_main"]
    custom_required: list = ["profile_files"]
    def __init__(self):
        super(HLM, self).__init__()
        self.data = {}

    def check_data(self, data: dict) -> None:
        super(HLM, self).check_data(data)

        hlm_type = data["hlm_type"]

        if hlm_type == "ramp-down":
            # Default ip_transition, no need to check
            pass
        elif hlm_type == "elm":
            super(HLM, self).check_data(data, self.elm_required)
        elif hlm_type == "custom":
            super(HLM, self).check_data(data, self.custom_required)
            # Now check if the files exist
            for file in data["profile_files"]:
                if not os.path.exists(file):
                    log.error(f"Plasma profile {file} does not exist!")
                    sys.exit(5)

        elif hlm_type == "single":
            super(HLM, self).check_data(data, self.single_required)
        elif hlm_type == "double":
            super(HLM, self).check_data(data, self.double_required)
        else:
            log.error(f"Wrong entry for hlm_type= {hlm_type}")
            sys.exit(3)

    def load_data(self, data: dict) -> None:
        self.check_data(data)

        self.set_required_attr(data)

        if self.data["hlm_type"] == "elm":
            self.set_required_attr(data, self.elm_required)
        elif self.data["hlm_type"] == "ramp-down":
            pass
        elif self.data["hlm_type"] == "single":
            self.set_required_attr(data, self.single_required)
        elif self.data["hlm_type"] == "double":
            self.set_required_attr(data, self.double_required)
        elif self.data["hlm_type"] == "custom":
            self.set_required_attr(data, self.custom_required)

        self.set_optional_attr(data)

class CASE(object):
    """Collection of the case data:
    """
    def __init__(self):
        import l2g.equil
        import l2g.comp

        # Input data
        self.geo_obj: GEOMETRY = None
        self.equ_obj: EQUILIBRIUM = None
        self.hlm_obj: HLM = None


        # Interface with L2G
        self.flt_obj: l2g.comp.FieldLineTracer = l2g.comp.FieldLineTracer()
        self.embree_shadow_geom_ids: list = []
        # For calcualting OMP objects.
        self.omp_obj: l2g.comp.FieldLineTracer = l2g.comp.FieldLineTracer()

        # For storing results
        self.med_obj: l2g.comp.MEDMeshIO = None
        # For collecting equilibria data
        self.ite_obj: l2g.equil.EquilibriumIterator = None

        # Additional variables
        # Output director
        self.output_directory: str = None
        self.result_file_path: str = None
        self.case_name: str = None
        self.run_flt: bool = False

        # Custom wall limiter points
        self.custom_wall_limiter: bool = False
        self.wall_limiter_r: list = False
        self.wall_limiter_z: list = False

        # Plasma displacement
        self.plasma_r_displ: float = 0.0
        self.plasma_z_displ: float = 0.0

        # Wall silhouette displacement
        self.wall_silh_r_shift: float = 0.0
        self.wall_silh_z_shift: float = 0.0

        # Align LCFS
        self.align_lcfs: bool = False
        # Custom LCFS values
        self.custom_lcfs_values: bool = False
        self.lcfs_values: List[float] = []

        # For obtaining FLs.
        self.fl_ids: list = None

        # Options
        self.get_field_lines: bool = True
        self.apply_heat_load: bool = True
        self.create_graphics: bool = True

        # Sub-graphics
        self.plot_heat_loads: bool = True
        self.plot_psi_maps: bool = True

    def set_objects(self, geo_obj: GEOMETRY, equ_obj: EQUILIBRIUM,
            hlm_obj: Optional[HLM]):

        self.geo_obj = geo_obj
        self.equ_obj = equ_obj
        self.hlm_obj = hlm_obj # Can be None

    def reset_flt_obj(self):
        import l2g.comp
        self.flt_obj: l2g.comp.FieldLineTracer = l2g.comp.FieldLineTracer()
        self.omp_obj: l2g.comp.FieldLineTracer = l2g.comp.FieldLineTracer()

    def load_target_mesh(self):
        """From the MEDMeshIO load the target data.
        """
        vertices, cells = self.med_obj.getMeshData()
        self.flt_obj.setTargetData(vertices, cells)

    def load_shadow_meshes(self):
        """From the
        """
        import l2g.utils.meshio
        import glob
        d = self.geo_obj.data
        flt = self.flt_obj

        shadowMeshFiles = []
        for fileName in d["shadow_meshes"]:
            if '*' in fileName:
                shadowMeshFiles += glob.glob(fileName)
            else:
                shadowMeshFiles.append(fileName)

        if "exclude_meshes" in d:
            mesh_to_remove = []
            # Since files are actual paths it is the easiest to just loop the list
            # and accumulate which meshes to remove
            for filePath in shadowMeshFiles:
                fileName = os.path.basename(filePath)
                if fileName in d["exclude_meshes"]:
                    mesh_to_remove.append(filePath)
            #
            for m in set(mesh_to_remove):
                shadowMeshFiles.remove(m)
        geom_ids = []

        if "include_target_in_shadow" in d:
            if d["include_target_in_shadow"]:
                if flt.target_vertices is not None:
                    log.info("Including target to shadowing Embree.")
                    geom_id = flt.embree_obj.commitMesh(
                        flt.target_vertices * flt.parameters.target_dim_mul,
                        flt.target_triangles)
                    geom_ids.append((geom_id, d["target_mesh"]))
                else:
                    log.info("There is no target data to include to Embree.")

        log.info(f"Loading {len(shadowMeshFiles)} mesh/es to Embree.")
        for filePath in shadowMeshFiles:
            fileName = os.path.basename(filePath)
            v, t = l2g.utils.meshio.readMesh(filePath)
            geom_id = flt.embree_obj.commitMesh(v * 1e-3, t)
            geom_ids.append((geom_id, filePath))

        if "afl_catcher_meshes" in d:
            afl_catcher_meshes = d["afl_catcher_meshes"]
            log.info(f"Loading {len(afl_catcher_meshes)} meshes, for catching " +
                     "and marking FLs as shadowed")
            for filePath in afl_catcher_meshes:
                v, t = l2g.utils.meshio.readMesh(filePath)
                geom_id = flt.embree_obj.commitMesh(v * 1e-3, t)
                flt.parameters.artificial_fl_catcher_geom_id.add(geom_id)
                fileName = os.path.basename(filePath)
                log.info(f"Loaded {fileName} as {geom_id}")

        log.info("Done.")
        self.embree_shadow_geom_ids = geom_ids

    def set_flt_objs(self):
        import l2g.workflow
        l2g.workflow.set_parameters_and_options(self.geo_obj.data,
                                                self.flt_obj)
        # Now set the mesh data.
        self.fl_ids = self.geo_obj.data["fl_ids"]
        # Propagate the HLM parameters to the fl object

        l2g.workflow.check_if_in_slurm_job(self.flt_obj)
        log.info("Parameters set to FieldLineTracer:")
        log.info(self.flt_obj.parameters.dump())
        log.info(self.flt_obj.options.dump())

        # Check align_lcfs
        if "align_lcfs" in self.geo_obj.data:
            flag = self.geo_obj.data["align_lcfs"]
            log.info(f"Align lcfs in case of limiter: {flag}")
            self.align_lcfs = flag

        if self.hlm_obj is None:
            self.apply_heat_load = False
            self.plot_heat_loads = False
            return

        self.flt_obj.hlm_params.load(self.hlm_obj.data)
        log.info(self.flt_obj.hlm_params.dump())

        if self.hlm_obj.data["hlm_type"] == "elm":
            self.omp_obj = l2g.comp.FieldLineTracer()
            l2g.workflow.load_elm_settings(self.hlm_obj.data, self.omp_obj)

            # Also set the same number of parallel threads.
            self.omp_obj.parameters.num_of_threads = \
                self.flt_obj.parameters.num_of_threads

            # Additionally check if run in SLURM job
            l2g.workflow.check_if_in_slurm_job(self.omp_obj)
            log.info("Parameters set for OMP FieldLineTracer")
            log.info(self.omp_obj.parameters.dump())
            log.info(self.omp_obj.options.dump())

    def create_target_medio(self):
        import l2g.comp
        self.med_obj = l2g.comp.MEDMeshIO()
        log.info(f"Reading target mesh data from: {self.geo_obj.data['target_mesh']}")
        self.med_obj.readMeshFromMedFile(self.geo_obj.data["target_mesh"])

    def prepare(self):
        # It just sets the input file
        self.create_target_medio()
        # From the target medio get the mesh data
        self.load_target_mesh()
        # Now load the shadow mesh data.
        self.load_shadow_meshes()
        # Set the options and parameters
        self.set_flt_objs()
        # Now the eqdsk
        self.prepare_equil_iterator_obj()

    def prepare_equil_iterator_obj(self):
        import l2g.equil
        self.ite_obj = l2g.equil.EquilibriumIterator()

        correct_helicity = True
        if "correct_helicity" in self.equ_obj.data:
            correct_helicity = self.equ_obj.data["correct_helicity"]
        self.ite_obj.correctHelicity(correct_helicity)

        if self.equ_obj.data["equilibrium_type"] == "eqdsk_files":
            self.ite_obj.loadEqdskEquilibriums(self.equ_obj.data["eqdsk_files"])
        elif self.equ_obj.data["equilibrium_type"] == "imas":
            print(self.equ_obj.data["imas"])
            self.ite_obj.loadIMASEquilibriums(self.equ_obj.data["imas"])

        if self.equ_obj.data["custom_wall_limiter"]:
            self.custom_wall_limiter = self.equ_obj.data["custom_wall_limiter"]
            log.info(f"Custom wall limiter: {self.custom_wall_limiter}")
            self.wall_limiter_r = self.equ_obj.data["wall_limiter_r"]
            self.wall_limiter_z = self.equ_obj.data["wall_limiter_z"]
            log.info(f"Loaded {len(self.wall_limiter_r)} R points")
            log.info(f"Loaded {len(self.wall_limiter_z)} Z points")

        if "plasma_r_displ" in self.equ_obj.data:
            self.plasma_r_displ = self.equ_obj.data["plasma_r_displ"]
        if "plasma_z_displ" in self.equ_obj.data:
            self.plasma_z_displ = self.equ_obj.data["plasma_z_displ"]

        if "wall_silh_r_shift" in self.equ_obj.data:
            self.wall_silh_r_shift = self.equ_obj.data["wall_silh_r_shift"]
        if "wall_silh_z_shift" in self.equ_obj.data:
            self.wall_silh_z_shift = self.equ_obj.data["wall_silh_z_shift"]
        self.ite_obj.applyWallSilhouetteShift(self.wall_silh_r_shift,
                                              self.wall_silh_z_shift)

    def see_if_results_exists(self, output_directory: str = "",
            med_file_name: str = ""):
        if output_directory == "":
            output_directory = os.getcwd()
        else:
            if os.path.exists(output_directory):
                # Check just to be sure if it is not a file
                if os.path.isfile(output_directory):
                    log.error(f"The output directory {output_directory} is a file!")
                    sys.exit(1)
            else:
                os.makedirs(output_directory)

        case_name = f'{self.geo_obj.data["name"]}_{self.equ_obj.data["name"]}'
        if med_file_name == "":
            # Overwrite the med_file_name AND case name so that all images
            # will consistently have the same base name.
            med_file_name = f"{case_name}.med"
        else:
            if not med_file_name.lower().endswith(".med"):
                case_name = med_file_name
                med_file_name = f'{med_file_name}.med'
            else:
                case_name = med_file_name[:-4]

        result_file_path = os.path.join(output_directory, med_file_name)

        # The additional files for connection length and elm profiles can have
        # indexes, regarding which equilibrium it is referencing.
        run_flt = False
        log.info(f"Checking {result_file_path}")
        if os.path.exists(result_file_path):
            log.info("FLT result file exists.")
            import medcoupling as mc
            field_names = mc.GetAllFieldNames(result_file_path)
            if "conlen" not in field_names:
                log.info("FLT result file exists but does not have conlen array.")
                run_flt = True
            else:
                iterations = mc.GetAllFieldIterations(result_file_path, "conlen")
                if len(iterations) != len(self.ite_obj):
                    run_flt = True
                    log.info("FLT result file does not have data for all " +
                                 "equilibriums. ")
        else:
            log.info("FLT result file doesn't exist.")
            run_flt = True

        # To be sure, let's also create other output files.

        log.info(f"Need to obtain FLT data?: {run_flt}")
        log.info(f"Name set for the case: {case_name}")
        log.info(f"Output directory set: {output_directory}")
        log.info(f"Name of result FLT MED file: {result_file_path}")

        self.run_flt = run_flt
        self.case_name = case_name
        self.result_file_path = result_file_path
        self.output_directory = output_directory

    def create_result_file_if_necessary(self):
        if self.run_flt:
            log.info(f"Writing target mesh to {self.result_file_path}")
            self.med_obj.writeMesh(self.result_file_path)
        else:
            # No FLT needed, but we still need to access the files.
            self.med_obj.med_file_path = self.result_file_path

    def main(self):
        """Main loop function for the case.
        """
        import l2g.comp

        import time
        time_start = time.perf_counter()
        N = len(self.ite_obj)

        hlm_type = self.flt_obj.hlm_params.hlm_type

        if self.apply_heat_load:
            # Print message which profile is applied
            if hlm_type == "elm":
                log.info("Applying Flat-Top plasma profile")
            elif hlm_type == "ramp-down":
                log.info("Applying Ramp-Down plasma profile")
            elif hlm_type == "single":
                log.info("Applying single exponential plasma profile")
            elif hlm_type == "double":
                log.info("Applying double exponential plasma profile")
            elif hlm_type == "custom":
                log.info("Applying custom profile")
            else:
                log.error(f"Unknown plasma profile: {hlm_type}. Stopping!")
                return
            import l2g.hlm.general
            import l2g.hlm.ramp_down
            import l2g.hlm.steady_state

        if self.create_graphics:
            # Additional arrays if they are needed
            time_array = []
            # Used for ramp-down graphics
            psol_array = []
            ip_array = []
            lambdaq_array = []
            # Create a simple X-axis drsep array for plotting Q_parallels.
            drsep = np.linspace(0, 0.3, 100)
            label = ""

            # Create path for file
            out_file = os.path.join(self.output_directory, self.case_name)

            hlm_type = self.flt_obj.hlm_params.hlm_type

            import matplotlib.pyplot as plt
            import matplotlib.figure
            import matplotlib.axes
            import l2g.plot

        for index, associated_time, equilibrium in self.ite_obj:
            log.info(f"Processing {index + 1} of {N} equilibriums. Time={associated_time:.3f}.")

            # Set the index and associated time to the med object.
            self.med_obj.setIndexAndTime(index, associated_time)

            # Adjust the equilibrium data if there is any custom points for the
            # wall silhouette.

            if self.custom_wall_limiter:
                equilibrium.wall_contour_r = np.array(self.wall_limiter_r)
                equilibrium.wall_contour_z = np.array(self.wall_limiter_z)

            # Apply the equilibrium to the FLT object
            self.flt_obj.setEquilibrium(equilibrium)
            # Apply the parameters so that the kernel gets activated and
            # properly configured.
            self.flt_obj.applyParameters()
            self.flt_obj.loadEq()

            # Paint the equilibrium quantities to the geometry
            self.flt_obj.processDataOnMesh()

            # If the equilibrium is limiter, then in case that the geometry is also
            # the one that has the LCFS contact point, we need to adjust r_move or
            # z_move in order to put the equilibrium data correctly on the geometry.

            # For now we shift only in the radial direction.
            self.flt_obj.evaluateEq()
            if self.custom_lcfs_values:
                # Assign the custom LCFS value
                self.flt_obj.eq.psiLCFS = self.lcfs_values[index]

            if self.flt_obj.eq.type_ == 'lim':
                if self.align_lcfs:
                    self.flt_obj.alignGeometryWithLCFS()

            if self.run_flt:
                # Run the FLT data since we have to run_flt.
                self.flt_obj.runFltOnMesh()
                # Calculate the drsep quantity on the geometry.
                self.flt_obj.calculateDrsep()
                # Now save the FLT data
                l2g.comp.dump_flt_mesh_results_to_med(
                    self.flt_obj.mesh_results, self.med_obj)

            else:
                # Now we already have the FLT results, so we need to load them
                l2g.comp.load_flt_mesh_results_from_med(
                    self.flt_obj.mesh_results, self.med_obj)

            # Optionals now.
            if self.get_field_lines:
                fl_base_name = ".".join(self.result_file_path.split(".")[:-1]) # without .med
                # Now get the FL
                if all(isinstance(x, list) for x in self.fl_ids):
                    # We have a group of IDs.
                    for j, subset in enumerate(self.fl_ids):
                        self.flt_obj.fl_ids = subset
                        self.flt_obj.getFL()
                        fl_path_name = os.path.join(self.output_directory,
                            f"{fl_base_name}_{j}_{index}.vtk")
                        l2g.comp.save_results_to_vtk(self.flt_obj.fl_results, fl_path_name)
                else:
                    self.flt_obj.fl_ids = self.fl_ids
                    self.flt_obj.getFL()
                    fl_path_name = os.path.join(self.output_directory,
                        f"{fl_base_name}_{index}.vtk")
                    l2g.comp.save_results_to_vtk(self.flt_obj.fl_results, fl_path_name)

            if self.apply_heat_load:
                if hlm_type == "elm":
                    self.omp_obj.setEquilibrium(equilibrium)
                    self.omp_obj.applyParameters()
                    self.omp_obj.loadEq()
                    self.omp_obj.evaluateEq()
                    self.obtain_elm_plm_profile(index)
                elif hlm_type == "custom":
                    self.load_custom_profile(index)

                self.flt_obj.applyHLM()

                log.info("Writing arrays to file")

                self.med_obj.writeArray(array=self.flt_obj.hlm_results.flux_expansion,
                    array_name="Total_flux expansion")
                self.med_obj.writeArray(array=self.flt_obj.hlm_results.q_inc,
                    array_name=r"$q_{inc}\;\;[\frac{W}{m^2}]$")
                self.med_obj.writeArray(array=self.flt_obj.hlm_results.q_par,
                    array_name=r"$q_{par}\;\;[\frac{W}{m^2}]$")

                if hlm_type == "elm":
                    # Writing additional arrays
                    # Take arrays from additional_arrays variable.
                    self.med_obj.writeArray(
                        array=self.flt_obj.hlm_results.additional_arrays[0],
                        array_name=r"$q_{\parallel}\;\;[\frac{W}{m^2}]$ ELM")
                    self.med_obj.writeArray(
                        array=self.flt_obj.hlm_results.additional_arrays[1],
                        array_name=r"$q_{\parallel}\;\;[\frac{W}{m^2}]$ inter-ELM")
                elif hlm_type == "ramp-down":
                    lambda_q = self.flt_obj.hlm_results.additional_arrays[0]
                    log.info(f'[RAMP DOWN]: i={index} t={associated_time:.3f} lq={lambda_q:.2f} Ip={self.flt_obj.equilibrium.Ip:.2e}')

            if self.create_graphics:
                time_array.append(associated_time)
                # Plot the psi
                base_name = os.path.join(self.output_directory,
                    f"{self.case_name}_{associated_time}")

                if self.plot_psi_maps:
                    figure = l2g.plot.plot_psi(eq=self.flt_obj.eq,
                                               output_path=base_name,
                                               save_files=True,
                                               plot_midplane=True)
                    # To avoid the "more than 20 figures open" matplotlib
                    # warning.
                    plt.close(figure)
                    del figure

                if self.plot_heat_loads:
                    # See if hlm is ELM.
                    if hlm_type == "elm":
                        self.omp_obj.setEquilibrium(equilibrium)
                        self.omp_obj.applyParameters()
                        self.omp_obj.loadEq()
                        self.omp_obj.evaluateEq()
                        self.obtain_elm_plm_profile(index)

                    if hlm_type == "custom":
                        self.load_custom_profile(index)

                    # Get IWL or OWL parameters
                    if self.flt_obj.parameters.side == "iwl":
                        Rb, _, Btotal, Bpm = self.flt_obj.eq.getIWL_midplane()
                    else: # owl
                        Rb, _, Btotal, Bpm = self.flt_obj.eq.getOWL_midplane()

                    if hlm_type == "custom":
                        label = "Custom"
                        q_par = l2g.hlm.general.custom(
                            drsep=drsep, points=self.flt_obj.hlm_params.points,
                            profile=self.flt_obj.hlm_params.profile)
                    elif  hlm_type == "elm":
                        label = "ELM PLM + inter ELM"
                        elm_par = l2g.hlm.general.custom(
                            drsep=drsep, points=self.flt_obj.hlm_params.points,
                            profile=self.flt_obj.hlm_params.profile)
                        interELM_par = l2g.hlm.steady_state.inter_ELM(
                            drsep=drsep, R_bdry=Rb, B_total=Btotal, B_pol=Bpm,
                            Rb=self.flt_obj.hlm_params.r_break)
                        q_par = elm_par + interELM_par
                    elif hlm_type == "single":
                        q_par = l2g.hlm.general.single_exponential_psol(
                            drsep=drsep, Bt=Btotal, Bpm=Bpm, Rb=Rb,
                            P_sol=self.flt_obj.hlm_params.p_sol,
                            F=0.5, lambda_q=self.flt_obj.hlm_params.lambda_q)
                        label = "Single exponential"
                    elif hlm_type == "double":
                        label = "Double exponential"
                        q_par = l2g.hlm.general.double_exponential_psol(
                            drsep=drsep, Bt=Btotal, Bpm=Bpm, Rb=Rb,
                            lambda_q_main=self.flt_obj.hlm_params.lambda_q_main,
                            lambda_q_near=self.flt_obj.hlm_params.lambda_q_near,
                            Rq=self.flt_obj.hlm_params.ratio,
                            P_sol=self.flt_obj.hlm_params.p_sol,
                            F=0.5)
                    elif hlm_type == "ramp-down":
                        label = "Ramp-Down"
                        psol_array.append(equilibrium.Psol)
                        ip_array.append(equilibrium.Ip)
                        # P_sol taken from IMAS, hopefully.
                        Ip = equilibrium.Ip
                        if Ip < self.flt_obj.hlm_params.ip_transition:
                            lambda_q = l2g.hlm.ramp_down.decay_length_L_mode_diverted(
                                a = equilibrium.a, Ip=Ip,
                                Area=equilibrium.Area, R=equilibrium.mag_axis_r)
                            q_par = l2g.hlm.general.single_exponential_psol(
                                drsep, Btotal, Bpm, Rb, lambda_q * 1e-3,
                                equilibrium.Psol)
                        else:
                            lambda_q = float("NaN")
                            q_par = np.zeros(drsep.shape)
                        lambdaq_array.append(lambda_q)



                    figure: matplotlib.figure.Figure = plt.figure()

                    ax: matplotlib.axes.Axes = figure.add_subplot()

                    # Plot the q_parallel against the drsep
                    ax.semilogy(drsep * 1e3, q_par, color="blue", label=label)
                    ax.grid(b=True, axis="both")
                    ax.set_xlabel(r"$r-r_{sep}[mm]$")
                    ax.set_title(self.case_name)

                    # Additional plots, in case of Steady-State
                    if hlm_type == "elm":
                        separatrix_distance = self.flt_obj.eq.drsep
                        ax.semilogy(drsep * 1e3, elm_par, 'y--', label="ELM PLM",)
                        ax.semilogy(drsep * 1e3, interELM_par, 'g--',
                            label="inter-ELM")
                        _drsep = np.linspace(separatrix_distance, 300, 100)
                        old_q_par = 5 * np.exp(-(_drsep - separatrix_distance) / 0.09)\
                                  + 3 * np.exp(-(_drsep - separatrix_distance) / 0.17)
                        ax.semilogy(_drsep, old_q_par, 'k')

                        ax.vlines(separatrix_distance, ymin=0, ymax=200, colors='r')
                        ax.text(separatrix_distance, 150, r"$2^{nd} sep$", rotation=90)

                    # Save the figure
                    figure.savefig(f"{out_file}_{index}_qpar.pdf")
                    figure.savefig(f"{out_file}_{index}_qpar.png")

                    figure.clf()
                    plt.close(figure)
                    del figure

        if self.create_graphics and self.plot_heat_loads and hlm_type == "ramp-down":
            psol_array = np.array(psol_array)
            ip_array = np.array(ip_array)
            # Plot the Ip, P_sol and lambda_q

            figure: matplotlib.figure.Figure = plt.figure()
            ax: matplotlib.axes.Axes = figure.add_subplot()
            ax_twin: matplotlib.axes.Axes = ax.twinx()

            # Plot on left Y-axis Psol and Ip
            p1 = ax.plot(time_array, psol_array / 1e6, color="g", label=r"P_{sol}")
            p2 = ax.plot(time_array, ip_array / 1e6, color="r", label=r"I_{p}")
            ax.set_ylabel(r"$I_p$ [MA], $P_{sol}$ [MW]")
            ax.set_xlabel("Time [s]")
            ax.grid()

            p3 = ax_twin.plot(time_array, lambdaq_array, color="b",
                              label=r"$\lambda_q$")
            plots = p1 + p2 + p3
            labels = [_.get_label() for _ in plots]
            ax.legend(plots, labels)
            figure.savefig(f"{out_file}_ip_psol_lambdaq.pdf")
            figure.savefig(f"{out_file}_ip_psol_lambdaq.png")
            figure.clf()
            plt.close(figure)
            del figure
        total_time = time.perf_counter() - time_start
        log.info(f"Finished in {total_time:.2f} seconds.")

    def obtain_elm_plm_profile(self, associated_time=""):
        """This function should only be called when calling applyHLM function,
        as the input equilibrium data has to be loaded to the ELM object.
        """
        # Pre-calculation step! To obtain ELM profile points and
        # profile
        # See if the ELM PLM profile file exist.

        # Generate the file path names
        import l2g.hlm.steady_state

        if associated_time != "":
            associated_time = f"_{associated_time}"

        elm_qpar_name = f"{self.case_name}{associated_time}_elm_qpar.dat"
        conlen_name = f"{self.case_name}{associated_time}_conlen.dat"
        elm_qpar_file_path = os.path.join(self.output_directory, elm_qpar_name)
        conlen_file_path = os.path.join(self.output_directory, conlen_name)

        if os.path.exists(elm_qpar_file_path):
            qelm_data = np.loadtxt(elm_qpar_file_path)
        else:
            # First check if connection length files exist.
            if os.path.exists(conlen_file_path):
                conlen_data = np.loadtxt(conlen_file_path)
            else:
                self.omp_obj.obtainOwlConlenGraph()
                conlen_data = self.omp_obj.owl_conlen_data
                np.savetxt(conlen_file_path, conlen_data,
                    header="drsep(Phi) [m], drsep(R) [m], Connection length " +
                           "Up [m], Connection length Down [m]")
            # Output name
            Rb, _, Btotal, Bpm = self.omp_obj.eq.getOWL_midplane()

            graphics_output_name = os.path.join(self.output_directory,
                f"{self.case_name}{associated_time}")
            qelm_data = l2g.hlm.steady_state.get_elm_data(conlen_data,
                generate_graphics=True, output_name=graphics_output_name,
                Rb=Rb, Btot=Btotal, Bpm=Bpm,
                r_break=self.omp_obj.parameters.r_break,
                drsep=self.omp_obj.eq.drsep)
            np.savetxt(elm_qpar_file_path, qelm_data,
                header="drsep [m], elm PLM [MW/m^2]")

        elm_data_r = qelm_data[:, 0]
        elm_data_q = qelm_data[:, 1] * 1e6

        # Now set the profile data to the main flt object.
        self.flt_obj.hlm_params.points = elm_data_r
        self.flt_obj.hlm_params.profile = elm_data_q

    def load_custom_profile(self, index):
        """Load a custom profile.
        """
        # See if we have a list of profile files.
        file = self.hlm_obj.data["profile_files"][index]
        data = np.loadtxt(file)
        if data.shape[0] == 2:
            self.flt_obj.hlm_params.points = data[0]
            self.flt_obj.hlm_params.profile = data[1]
        else:
            self.flt_obj.hlm_params.points = data[:, 0]
            self.flt_obj.hlm_params.profile = data[:, 1]
