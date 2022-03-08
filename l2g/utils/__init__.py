# Main function that reads mesh data
# from ._meshio import readMesh
# MEDCoupling function for transforming numpy array to field and vice versa
# from ._meshio import fieldToNumpy, numpyArrayToField, writeFieldToAlreadyExistingMesh

import l2g.utils.meshio

# Next functions are helper functions for
import glob
import re
import json
import logging
import os
import sys
log = logging.getLogger(__name__)

import numpy as np


import typing

if typing.TYPE_CHECKING:
    from l2g.comp import FieldLineTracer

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
        [log.info(line) for line in lines[e.lineno - 2: e.lineno]]
        log.info(lines[e.lineno] + '    <----- ERROR')
        [log.info(line) for line in lines[e.lineno + 1: e.lineno + 3]]
        raise

    return out



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


def load_l2g_json(file_path: str) -> dict:
    """Normal JSON format does not allow // comments, so we first clean it of
    such comment blocks and also remove empty lines.

    If there is a problem with the JSON file, lines are printed to show user,
    what is causing the error. Also since it fails to load the JSON file, do
    not allow the program to run.
    """

    # Load json file
    inp = {}
    text = open(file_path, 'r').read()
    inp = json_loads(text)
    return inp

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

    set_parameters_and_options(d, flt)

    if not os.access(d["target_mesh"], os.R_OK):
        log.error(f"Failed to read target mesh {d['target_mesh']}!")
        log.info("Check the validity of the path.")
        sys.exit(-1)

    # Load target meshes
    log.info(f"Loading target mesh data {d['target_mesh']}")
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
    geomIds = []
    if "include_target_in_shadow" in d:
        log.info("Including target to shadowing Embree.")
        if d["include_target_in_shadow"]:
            geomId = flt.embree_obj.commitMesh(
                verticesTarget * flt.parameters.target_dim_mul,
                trianglesTarget)
            geomIds.append((geomId, d["target_mesh"]))

    log.info(f"Loading {len(shadowMeshFiles)} mesh/es to Embree.")
    for filePath in shadowMeshFiles:
        fileName = os.path.basename(filePath)
        v, t = l2g.utils.meshio.readMesh(filePath)
        geomId = flt.embree_obj.commitMesh(v * 1e-3, t)
        geomIds.append((geomId, filePath))
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
