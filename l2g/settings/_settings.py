from typing import Hashable, Any
import logging

log = logging.getLogger(__name__)

class BaseClassOptions:
    __slots__: list[str] = []
    def dump(self) -> dict[Hashable, Any]:
        """Dumps a dictionary of the attributes defined in __slots__.

        Returns:
            out (dict): Dictionary of all the attributes and the their values
                defined in __slots__.
        """
        out: dict[Hashable, Any] = {}
        for slot in self.__slots__:
            out[slot] = self.__getattribute__(slot)
        return out

    def load(self, options: dict[str, Any]) -> None:
        """Loads values from input dictionary to __slots__.

        Argument:
            options (dict): Dictionary of attributes and their values. If an
                attribute is not defined

        """
        for option in options:
            if option in self.__slots__:
                self.__setattr__(option, options[option])
            else:
                log.warning(f"Option {option} not defined in {self.__class__.__name__}. Skipping")
        return

class Options(BaseClassOptions):
    """Options on how to run FLT. These settings are implemention wise, meaning
    they determine how the underlying kernel is called.
    """
    __slots__: list[str] = ["switch_getFL_with_FLT", "switch_runFLT",
                            "switch_applyPP"]

    def __init__(self):
        #: getFL_switch:
        #: False - get FL with FLT off
        #: True - get FL with FLT on
        self.switch_getFL_with_FLT: bool = True

        #: switch_runFLT:
        #: 0 - trace until end of some toroidal angle (t_end)
        #: 1 - trace until some length (max_conlen).
        self.switch_runFLT: int = 1

        #: switch_applyPP:
        #: 0 - No plasma profile
        #: 1 - Single Exponential profile
        self.switch_applyPP: int = 0

class Parameters(BaseClassOptions):
    __slots__: list[str] = ["plasma_r_displ", "plasma_z_displ", "wall_silh_r_displ",
                            "wall_silh_z_displ", "time_end", "time_step",
                            "max_fieldline_length",
                            "self_intersection_avoidance_length", "abs_error",
                            "rel_error", "target_to_m", "shadow_to_m",
                            "num_of_threads", "side", "cutoff_conlen",
                            "artificial_fl_catcher_geom_id",
               "lcfs_max_align_dist"]
    def __init__(self):

        #: The plasma_*_displ parameters displaces the plasma in either R or Z
        #: direction. In meters.
        self.plasma_r_displ: float = 0.0
        self.plasma_z_displ: float = 0.0

        #: Equilibrium wall silhouette displacement. In meters. Use when the
        #: wall does not align with the used geometry.
        self.wall_silh_r_displ: float = 0.0
        self.wall_silh_z_displ: float = 0.0

        #: Time stop is the parameter for setting the toroidal angle end in
        #: radians.
        self.time_end: float = 1.0
        #: Time step is the resolution at which the solver checks the FL path
        #: for intersections. In other words toroidal resolution for solver.
        #: In radians.
        self.time_step: float = 0.01
        #: Maximum connection length when activating traces where we want to
        #: follow a FL until the desired length. In meters.
        self.max_fieldline_length: float = 100.0

        #: Self intersection avoidance length is a length in which the solver
        #: follows a FL slowly for a brief distance and then actually starts the
        #: FLT from the new origin point. In meters
        self.self_intersection_avoidance_length: float = 0.001

        #: Desired accuracy
        self.abs_error: float = 1e-4
        self.rel_error: float = 1e-4

        #: The following target_to_m is used to transform the units of mesh
        #: data to meters.
        self.target_to_m: float = 1e-3

        #: The following shadoDimMul is used to transform the units of
        #: coordinates of a ray when ray tracing into the same units that is the
        #: shadow.
        #: I.e. : shadowing geometry has mm, while rays are in meters, therefore
        #: we multiple the points of a ray to meters by multiplying with this
        #: factor.
        self.shadow_to_m: float = 1e-3

        #: Parameter for specifying number of CPU threads. For use with Cython
        #: OpenMP.
        self.num_of_threads: int = 1

        #: Parameter side for determining for which midplane, inner or outer.
        #: For use with HLM. Either "iwl" or "owl"
        self.side: str = "iwl"

        #: Advanced parameter for aligning LCFS to input 3D geometry. This is
        #: used when large encompassing 3d geometries are used with small
        #: limiter plasmas. In meters.
        self.lcfs_max_align_dist: float | int = 1

        #: Parameter cutoff_conlen specifies the length at which we consider
        #: that a field line wetts target area. In meters.
        self.cutoff_conlen: float | int = 4.3


        #: IDs of geometries that marks fieldlines as shadowed if they
        #: intercept a field line. Basically using geometry that encompasses
        #: the plasma volume and is used as a boundary where there is no plasma
        #: (example, vacuum vessel), then this array of IDs is used when
        #: applying heat load maps or evaluating the wetted area. If the ID
        #: of intercepted geometry belongs to this array for a given fieldline,
        #: then mark the fieldline as non-wetting on the target.
        self.artificial_fl_catcher_geom_id: set[int] = set()

class HLM(BaseClassOptions):
    """A common class that holds the parameters used in the implemented
    heat loads.
    """

    # Ramp - Down
    __slots__ = ["hlm_type", "p_sol", "lambda_q", "lambda_q_main",
                 "lambda_q_near", "ratio", "r_break", "ip_transition",
                 "Rb", "Z", "Btotal", "points", "profile", "additional_points",
                 "additional_profiles", "extrapolate", "outside_value",
                 "longwave_misaligment_applied", "longwave_l", 'q_par0']

    # Lists for each scenario
    flat_top: list[str] = ["r_break", "p_sol", "lambda_q_main",
                           "lambda_q_near"]
    ramp_down: list[str] = ["ip_transition"]
    single_exp: list[str] = ["p_sol", "lambda_q"]
    double_exp: list[str] = ["p_sol", "lambda_q_main", "lambda_q_near",
                             "ratio"]
    custom: list[str] = ["points", "profile", "additional_points",
                         "additional_profiles", "extrapolate", "outside_value",
                         "longwave_misaligment_applied", "longwave_l"]

    acceptable_types: list[str] = ["flat_top", "ramp_down", "single_exp",
                                   "double_exp", "custom"]
    # L-mod flat-top is similar to H-mod flat-top, except we use the scaling
    # law to determine the lambda_q_near and also use it as a r_break!
    L_mode: list[str] = ["r_break", "p_sol", "lambda_q_near", "lambda_q_main"]

    def __init__(self):
        self.hlm_type: str = ""
        self.p_sol: float = 100e6
        self.q_par0: float = 4e6
        self.lambda_q: float = 0.012
        self.lambda_q_main: float = 0.17
        self.lambda_q_near: float = 0.005
        self.r_break: float = 0.025
        self.ratio: float = 4
        self.ip_transition: float = 10e6

        self.points: list[float] = []
        self.profile: list[float] = []
        self.additional_points: list[float] = []
        self.additional_profiles: list[float] = []
        self.extrapolate: bool = True
        self.outside_value: float = 0.0
        self.longwave_misaligment_applied: bool = False
        self.longwave_l: float = 0.0
        # Not user set
        self.Rb: float = 0
        self.Z: float = 0 #
        self.Btotal: float = 0.0


if __name__ == "__main__":

    x = Parameters()
    d = x.dump()

    for key in d:
        print(key)