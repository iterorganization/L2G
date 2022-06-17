from l2g.workflow._case import GEOMETRY, EQUILIBRIUM, HLM, CASE

import yaml
import re
import json

from typing import List, Tuple, Optional, TYPE_CHECKING
import logging
log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from l2g.comp import FieldLineTracer

def load_yaml_configuration(yaml_text: str) -> Tuple[List[GEOMETRY], List[EQUILIBRIUM], List[HLM]]:
    """Reads a YAML configuration of a collection of input data for FLT cases
    and returns 3 lists.

    Arguments:
        yaml_text (str): String containing YAML documents.

    Returns:
        tuple: A tuple of 3 lists in the following order. A list of geometries
            or in other words the input geometry data and parameters for the
            FLT component, a list of equilibriums and a list of heat load
            mappings.
    """

    data = yaml.safe_load_all(yaml_text)

    geometry_objects: List[GEOMETRY] = []
    equilibrium_objects: List[EQUILIBRIUM] = []
    hlm_objects: List[HLM] = []

    for block in data:
        # Do necessary checks like type
        if "type" not in block:
            log.error("Block in the YAML has no type attribute, therefore skipping")
            log.debug(f"Block data: {block}")
            continue

        block_type = block["type"]
        # print(block)
        # Now create the objects and store them.
        if block_type == "geometry":
            geometry_obj = GEOMETRY()
            geometry_obj.load_data(block)
            geometry_objects.append(geometry_obj)
        elif block_type == "equilibrium":
            equilibrium_obj = EQUILIBRIUM()
            equilibrium_obj.load_data(block)
            equilibrium_objects.append(equilibrium_obj)
        elif block_type == "hlm":
            hlm_obj = HLM()
            hlm_obj.load_data(block)
            hlm_objects.append(hlm_obj)

    return geometry_objects, equilibrium_objects, hlm_objects

def json_remove_comments(json_like_string: str) -> str:
    r"""Removes C-style comments from *json_like* and returns the result.

    Regex pattern:

    ``//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"``

     - ``//.*?`` Match the comment section non-greedy zero or more times

    or

     - ``/\*.*?\*/`` Match the /* ... */ non-greedy zero or more times.

    The following matches find strings that contain also supposedly C-style
    comments. These are ignored then in the replacing function, as they are
    valid JSON strings.

     -  ``\'(?:\\.|[^\\\'])*\'`` Find strings that contain C-style comments.
     - ``"(?:\\.|[^\\"])*"`` Find strings that contain C-style comments.


    Example:

    .. code-block:: python

       test_json = '''\
       {
           "foo": "bar", // This is a single-line comment
           "baz": "blah" /* Multi-line
           Comment */
       }'''
       remove_comments('{"foo":"bar","baz":"blah",}')
         '{\n    "foo":"bar",\n    "baz":"blah"\n}'
    """
    comments_re = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    def replacer(match):
        s = match.group(0)
        if s[0] == '/': return ""
        return s
    return comments_re.sub(replacer, json_like_string)

def json_remove_trailing_commas(json_line_string: str) -> str:
    """Removes trailing commas from *json_like* and returns the result.

    Example:

    .. code-block:: python

       remove_trailing_commas('{"foo":"bar","baz":["blah",],}')
       '{"foo":"bar","baz":["blah"]}'
    """
    trailing_object_commas_re = re.compile(
        r'(,)\s*}(?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    trailing_array_commas_re = re.compile(
        r'(,)\s*\](?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    # Fix objects {} first
    objects_fixed = trailing_object_commas_re.sub("}", json_line_string)
    # Now fix arrays/lists [] and return the result
    return trailing_array_commas_re.sub("]", objects_fixed)

def json_loads(json_like: str) -> dict:
    """JSON format does not allow comments or trailing commas, therefore this
    functions first cleans the JSON like content and finally calls json.loads.
    """

    # Remove comments
    almost_json = json_remove_comments(json_like)

    # Remove trailing white spaces
    proper_json = json_remove_trailing_commas(almost_json)

    out = None

    try:
        out = json.loads(proper_json)
    except json.JSONDecodeError as e:
        # Get the line number error
        lines = proper_json.splitlines()
        log.info("Error when parsing JSON.")
        [log.info(line) for line in lines[e.lineno - 2: e.lineno-1]]
        log.info(lines[e.lineno-1] + '    <----- ERROR')
        [log.info(line) for line in lines[e.lineno: e.lineno + 3]]
        raise

    return out

def load_json_configuration(json_text: str) -> Tuple[List[GEOMETRY], List[EQUILIBRIUM], List[HLM]]:
    # First get the JSON structure
    data: dict = json_loads(json_text)

    geometry_objects: List[GEOMETRY] = []
    equilibrium_objects: List[EQUILIBRIUM] = []
    hlm_objects: List[HLM] = []
    # Now parse the data

    case_name = data["name"]

    # Load the FLT block
    geo_obj = GEOMETRY()
    geo_dict: dict = {}
    if "cutoff" in data:
        geo_dict["cutoff"] = data["cutoff"]
    if "fl_ids" in data:
        geo_dict["fl_ids"] = data["fl_ids"]

    if "parameters" in data:
        geo_dict["parameters"] = data["parameters"]
    if "options" in data:
        geo_dict["options"] = data["options"]

    geo_dict.update(data["flt"])
    geo_dict["name"] = case_name

    geo_obj.load_data(geo_dict)
    geometry_objects.append(geo_obj)

    # Load the equilibrium data
    equ_obj = EQUILIBRIUM()

    equ_dict = {}
    equ_dict["name"] = case_name
    if data["eq_type"] == "eqdsk":
        equ_dict["equilibrium_type"] = "eqdsk_files"
        equ_dict["eqdsk_files"] = data["eqdsk_files"]
    else:
        equ_dict["equilibrium_type"] = "imas"
        equ_dict["imas"] = data["imas"]
    if "wall_limiter" in data:
        equ_dict["custom_wall_limiter"] = True
        equ_dict["wall_limiter_r"] = data["wall_limiter"]["r"]
        equ_dict["wall_limiter_z"] = data["wall_limiter"]["z"]

    equ_obj.load_data(equ_dict)
    equilibrium_objects.append(equ_obj)
    # Load any or all HLM
    if "elm" in data:
        hlm_obj = HLM()
        hlm_dict: dict = {}
        hlm_dict["name"] = "elm"
        hlm_dict["hlm_type"] = "elm"
        hlm_dict.update(data["elm"])
        hlm_obj.load_data(hlm_dict)
        hlm_objects.append(hlm_obj)

    if "ramp-down" in data:
        hlm_obj = HLM()
        hlm_dict: dict = {}
        hlm_dict["name"] = "ramp-down"
        hlm_dict["hlm_type"] = "ramp-down"
        if "Ip transition" in data["ramp-down"]:
            hlm_dict["ip_transition"] = data["ramp-down"]["Ip transition"]
        hlm_obj.load_data(hlm_dict)
        hlm_objects.append(hlm_obj)
    if "single" in data:
        hlm_obj = HLM()
        hlm_dict: dict = {}
        hlm_dict["name"] = "single"
        hlm_dict["hlm_type"] = "single"
        hlm_dict[""]
        hlm_obj.load_data(hlm_dict)
        hlm_objects.append(hlm_obj)
    if "double" in data:
        hlm_obj = HLM()
        hlm_dict: dict = {}
        hlm_dict["name"] = "double"
        hlm_dict["hlm_type"] = "double"
        hlm_obj.load_data(hlm_dict)
        hlm_objects.append(hlm_obj)
    for key in data:
        if key.startswith("custom-hlm"):
            hlm_obj = HLM()
            hlm_dict: dict = {}
            hlm_dict["name"] = key[10:].replace("-", "_").lower()
            hlm_dict["hlm_type"] = "custom"
            hlm_dict.update(data[key])
            hlm_obj.load_data(hlm_dict)
            hlm_objects.append(hlm_obj)

    return geometry_objects, equilibrium_objects, hlm_objects

def create_case(geo_obj: GEOMETRY, equ_obj: EQUILIBRIUM,
        hlm_obj: Optional[HLM]) -> CASE:
    """Function that creates a case obj and ties the input data objects to it.

    Arguments:
        geo_obj (GEOMETRY): Contains information about geometry input data for
            the FLT case
        equ_obj (EQUILIBRIUM): Contains information about equilibrium input
            data
        hlm_obj (Optional[HLM]): Contains information about the heat load map

    Returns:
        case_obj (CASE): Object containing pre-defined functions for running a
            L2G flt case and its other utility functions.
    """
    case_obj = CASE()

    case_obj.geo_obj = geo_obj
    case_obj.equ_obj = equ_obj
    if hlm_obj is not None:
        case_obj.hlm_obj = hlm_obj
    return case_obj

def check_if_in_slurm_job(*flt_objs):
    """This function checks if the process is run in a SLURM environment. A
    simple check if environment variable SLURM_JOBID exists. If it exists, then
    it is deduced that the process is run inside a compute node, therefore
    ideally we increase the number of max OpenMP threads to whatever is the
    maximum for the compute node. One thing less to worry about the user.


    Arguments:
        flt_objs (FieldLinesTracer): Set of FLT objects, for which we set the
            maximum number of cpu threads.
    """
    import os
    if "SLURM_JOBID" in os.environ:
        log.info('Detected SLURM environment.')
        job_id = os.environ['SLURM_JOBID']
        log.info(f"SLURM JOB ID: {job_id}")
        log.info("Setting the number of OpenMP threads equal to the number of "
              + "allocated CPUs")
        # In this case assign the number of threads equal to the number of
        # CPUs per task.
        # We only have 1 distributed task.

        if "SLURM_CPUS_PER_TASK" in os.environ:
            log.info('Trying to get maximum number of cpus assigned for this task.')
            try:
                cpus_per_task = int(os.environ["SLURM_CPUS_PER_TASK"])
                for flt_obj in flt_objs:
                    flt_obj.parameters.num_of_threads = cpus_per_task
                log.info(f"case.parameters.num_of_threads={cpus_per_task}")
            except:
                log.info("Failed to obtain a number from SLURM_CPUS_PER_TASK")

def set_parameters_and_options(d: dict, flt_obj) -> None:
    """Set parameters to a FieldLineTracer object from a dictionary.
    """
    # Set the parameters. If unknown, just print and ignore
    if "parameters" in d:
        for parameter in d['parameters']:
            if not hasattr(flt_obj.parameters, parameter):
                log.info(f"Illegal parameter: {parameter}. Ignored")
            else:
                setattr(flt_obj.parameters, parameter, d['parameters'][parameter])

    # Set the options
    if "options" in d:
        for option in d['options']:
            if not hasattr(flt_obj.options, option):
                log.info(f"Illegal option: {option}. Ignored")
            else:
                setattr(flt_obj.options, option, d['options'][option])


def load_flt_settings(d: dict, flt) -> None:
    """Loads the settings for the FLT from a dictionary. Also loads geometries.


    .. todo::

       Figure out whether it is okay to have a function like this which loads
       settings and ALSO loads geometry.

    """
    import os
    import sys
    import glob

    set_parameters_and_options(d, flt)

    if not os.access(d["target_mesh"], os.R_OK):
        log.error(f"Failed to read target mesh {d['target_mesh']}!")
        log.info("Check the validity of the path.")
        sys.exit(-1)

    # Load target meshes
    log.info(f"Loading target mesh data {d['target_mesh']}")
    import l2g.utils.meshio
    verticesTarget, trianglesTarget = l2g.utils.meshio.readMesh(d["target_mesh"])

    flt.setTargetData(verticesTarget, trianglesTarget)
    log.info("Loading shadow mesh data.")
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
        log.info("Including target to shadowing Embree.")
        if d["include_target_in_shadow"]:
            geom_id = flt.embree_obj.commitMesh(
                verticesTarget * flt.parameters.target_dim_mul,
                trianglesTarget)
            geom_ids.append((geom_id, d["target_mesh"]))

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

def load_elm_settings(d: dict, flt: 'FieldLineTracer') -> None:
    """Prepares a FieldLineTracer object for evaluating ELM profiles.

    Essentially, in this case the main goal is to obtain the outer wall
    midplane connection length graph. For that what we really need is that the
    maximum length until which a field line is followed is set to a high value
    and that's it. This is then used for evaluating the ELM contribution via
    the ELM PLM model.

    """
    flt.parameters.time_step = 0.01
    flt.parameters.time_end = 2 * 3.141592653
    flt.parameters.max_connection_length = 1000
    flt.options.switch_runFLT = 1
    flt.parameters.target_dim_mul = 1

    set_parameters_and_options(d, flt)

    # embreeObj = l2g.comp.core.PyEmbreeAccell()

    shadowMeshFiles = d['shadow_meshes']
    flt.commitMeshesToEmbree(shadowMeshFiles)


    # for f in shadowMeshFiles:
    #     v, t = l2g.utils.meshio.readMesh(f)
    #     embreeObj.commitMesh(v * 1e-3, t)
    #     flt.setEmbreeObj(embreeObj)
