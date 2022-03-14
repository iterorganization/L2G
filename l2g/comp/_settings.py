"""These are not global settings for the whole package. This file contains
classes for setting a FLT study.
"""
import logging
import os
from pathlib import Path
from l2g.utils import load_l2g_json

log = logging.getLogger(__name__)

class Options:
    """Options on how to run FLT. These settings are implemention wise, meaning
    they determine how the underlying kernel is called.
    """
    __slots__=["switch_getFL_with_FLT", "switch_runFLT", "switch_applyPP"]
    def __init__(self):
        #: getFL_switch:
        #: False - get FL with FLT off
        #: True - get FL with FLT on
        self.switch_getFL_with_FLT = True

        #: switch_runFLT:
        #: 0 - trace until end of some toroidal angle (t_end)
        #: 1 - trace until some length (max_conlen).
        self.switch_runFLT = 0

        #: switch_applyPP:
        #: 0 - No plasma profile
        #: 1 - Single Exponential profile
        self.switch_applyPP = 0


class Parameters:
    __slots__=["plasma_r_displ", "plasma_z_displ", "wall_silh_r_displ",
               "wall_silh_z_displ", "time_end",
               "time_step", "max_connection_length",
               "self_intersection_avoidance_length", "abs_error", "rel_error",
               "target_dim_mul", "shadow_dim_mul", "num_of_threads",
               "side", "P_sol", "F_split", "q_parallel", "lambda_q_main",
               "lambda_q_near", "R_q", "rd_ip_transition",
               "cutoff_conlen", "r_break", "artificial_fl_catcher_geom_id"]
    def __init__(self):

        #: The plasma_*_displ parameters displaces the plasma in either R or Z
        #: direction. In meters.
        self.plasma_r_displ = 0.0
        self.plasma_z_displ = 0.0

        #: Equilibrium wall silhouette displacement. In meters. Use when the
        #: wall does not align with the used geometry.
        self.wall_silh_r_displ = 0.0
        self.wall_silh_z_displ = 0.0

        #: Time stop is the parameter for setting the toroidal angle end in
        #: radians.
        self.time_end = 1.0
        #: Time step is the resolution at which the solver checks the FL path
        #: for intersections. In other words toroidal resolution for solver
        self.time_step = 0.01
        #: Maximum connection length when activating traces where we want to
        #: follow a FL until the desired length. In meters.
        self.max_connection_length = 1.0

        #: Self intersection avoidance length is a length in which the solver
        #: follows a FL slowly for a brief distance and then actually starts the
        #: FLT from the new origin point. In meters
        self.self_intersection_avoidance_length = 0.005

        #: Desired accuracy
        self.abs_error = 1e-4
        self.rel_error = 1e-4

        #: The following targetDimMul is used to transform the units of mesh
        #: data to meters.
        self.target_dim_mul = 1e-3

        #: The following shadoDimMul is used to transform the units of
        #: coordinates of a ray when ray tracing into the same units that is the
        #: shadow.
        #: I.e. : shadowing geometry has mm, while rays are in meters, therefore
        #: we multiple the points of a ray to meters by multiplying with this
        #: factor.
        self.shadow_dim_mul = 1000

        #: Parameter for specifying number of CPU threads. For use with Cython
        #: OpenMP.
        self.num_of_threads = 1

        #: Parameter side for determining for which midplane, inner or outer.
        #: For use with HLM. Either "iwl" or "owl"
        self.side = "iwl"

        #: Power @SOL parameter for use with HLM
        self.P_sol = 1.0

        #: Split parameter F for use with HLM.
        self.F_split = 0.5

        #: Parallel value of heat load at midplane for use with HLM
        self.q_parallel = None

        #: Main decay length parameter for use with HLM
        self.lambda_q_main = 0.012

        #: Near decay length parameter for use with HLM.
        self.lambda_q_near = 0.005

        #: Double exponential parameter ratio Rq
        self.R_q = 1

        #: Double exponential parameter breakpoint r_break. In meters
        self.r_break = 0.025

        #: Ramp down transition
        self.rd_ip_transition = 10e6

        #: Parameter cutoff_conlen specifies the length at which we consider
        #: that a field line wetts target area. In meters.
        self.cutoff_conlen = 4.3e3


        #: IDs of geometries that marks fieldlines as shadowed if they
        #: intercept a field line. Basically using geometry that encompasses
        #: the plasma volume and is used as a boundary where there is no plasma
        #: (example, vacuum vessel), then this array of IDs is used when
        #: applying heat load maps or evaluating the wetted area. If the ID
        #: of intercepted geometry belongs to this array for a given fieldline,
        #: then mark the fieldline as non-wetting on the target.
        self.artificial_fl_catcher_geom_id = set()