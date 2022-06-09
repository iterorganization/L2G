from l2g.workflow._case import GEOMETRY, EQUILIBRIUM, HLM, CASE

import yaml
import re
import json

from typing import List, Tuple, Optional
import logging
log = logging.getLogger(__name__)

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