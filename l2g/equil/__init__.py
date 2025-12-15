from ._equilibrium import Equilibrium, correct_equilibrium_helicity
from ._eqdskg import EQDSKIO
from l2g.external.equilibrium_analysis import EQA as EQ
from ._owl_connection_graph import getOwlConlensGraph

from typing import Any

def getEquilibriumFromEQFile(eqdsk_file: str,
                             correct_helicity: bool=True) -> Equilibrium:
    """Gather the data required from an EQDSKIO object and populate the
    Equilibrium class.
    """

    eqdsk = EQDSKIO(eqdsk_file)
    return getEquilibriumFromEQDSKG(eqdsk, correct_helicity)

def getEquilibriumFromEQDSKG(eqdsk_obj: EQDSKIO,
                             correct_helicity: bool=True) -> Equilibrium:
    """Gather the data required from an EQDSKIO object and populate the
    Equilibrium class.
    """
    import numpy as np
    obj = Equilibrium()

    # Write the wall Limiter
    obj.wall_contour_r = eqdsk_obj.getRLIM()
    obj.wall_contour_z = eqdsk_obj.getZLIM()

    # Write the magnetic axis position
    obj.mag_axis_r = eqdsk_obj.RMAXIS
    obj.mag_axis_z = eqdsk_obj.ZMAXIS


    # Write the R-Z grid
    RLEFT, ZMID = eqdsk_obj.RLEFT, eqdsk_obj.ZMID
    RDIM, ZDIM = eqdsk_obj.RDIM, eqdsk_obj.ZDIM
    NW, NH = eqdsk_obj.NW, eqdsk_obj.NH
    ZMIN = ZMID - 0.5 * ZDIM
    ZMAX = ZMID + 0.5 * ZDIM

    obj.grid_r = np.linspace(RLEFT, RLEFT + RDIM, NW)
    obj.grid_dim_r = NW

    obj.grid_z = np.linspace(ZMIN, ZMAX, NH)
    obj.grid_dim_z = NH

    # Poloidal magnetic flux
    obj.psi = np.asarray(eqdsk_obj.PSIRZ)
    obj.psi_boundary = eqdsk_obj.SIBRY
    obj.psi_axis = eqdsk_obj.SIMAG
    obj.Ip = eqdsk_obj.CURRENT

    obj.psi_sign = np.sign(obj.psi_boundary - obj.psi_axis )

    # Write the Plasma current
    # Write the FPOL
    obj.fpol = np.asarray(eqdsk_obj.FPOL)
    obj.fpol_flux = np.linspace(eqdsk_obj.SIMAG, eqdsk_obj.SIBRY, NW)
    # Write the FPOL in vacuum. Just take the last value
    obj.fpol_vacuum = obj.fpol[-1]

    # Equilibrium type
    obj.type = ''

    if correct_helicity:
        correct_equilibrium_helicity(obj)

    return obj

def getBackupIMASWallIds(shot=116000, run=4, database="ITER_MD",
                         user_name="public", backend="mdsplus",
                         version='3') -> Any:
    """Gets the default Wall IDS machine description from the ITER ITER_MD
    machine description database. 116000/4
    """

    import imas
    uri = f"imas:{backend}?user={user_name};{shot=};{run=};database={db_name};version={version}"
    db_entry = imas.DBEntry(uri, 'r')
    return db_entry.get("wall", autoconvert=False)

def getEquilibriumFromIMASSlice(shot: int, run: int, user: str = "public",
        machine: str = "iter", time: float = 0.0, version: str = "3",
        backend: str = "mdsplus", correct_helicity: bool = True) -> Equilibrium:

    import imas
    import imas.ids_defs

    uri = f"backend:{backend}?user={user};{shot=};{run=};database={machine};version={version}"
    db_entry = imas.DBEntry(uri, 'r')
    summary_ids = db_entry.get("summary", autoconvert=False)
    wall_ids = db_entry.get("wall", autoconvert=False)

    equilibrium_ids = db_entry.get_slice("equilibrium", time,
                                         imas.ids_defs.CLOSEST_INTERP,
                                         autoconvert=False)

    equilibrium =  getEquilibriumFromIMAS(equilibrium_ids.time_slice[0],
        equilibrium_ids.vacuum_toroidal_field,
        wall_ids, summary_ids, correct_helicty=correct_helicity)
    return equilibrium

def getEquilibriumFromIMAS(equilibrium_ids_time_slice: Any,
                           vacuum_toroidal_field_ids: Any, wall_ids: Any,
                           summary_ids: Any=None,
                           correct_helicty: bool=True) -> Equilibrium:
    """Populate the Equilibrium class with relevant data.

    It's best to work on time slices of data, that means that
    equilibrium.time_slice is a single time slice, same for
    vacuum_toroidal_field.

    Optional:
        summary_ids can be optionally provided to obtain other relevant values.
    """
    import numpy as np
    obj = Equilibrium()
    try:
        # Write the wall limiter
        obj.wall_contour_r = list(wall_ids.description_2d[0].limiter.unit[0].outline.r)
        obj.wall_contour_z = list(wall_ids.description_2d[0].limiter.unit[0].outline.z)
    except:
        # Soooooo let's try the backup wall...
        wall_backup = getBackupIMASWallIds()

        wall_r = np.concatenate([wall_backup.description_2d[0].limiter.unit[0].outline.r,
                                 wall_backup.description_2d[0].limiter.unit[1].outline.r[::-1]])
        wall_z = np.concatenate([wall_backup.description_2d[0].limiter.unit[0].outline.z,
                                 wall_backup.description_2d[0].limiter.unit[1].outline.z[::-1]])
        obj.wall_contour_r = list(wall_r)
        obj.wall_contour_z = list(wall_z)

    # Write the magnetic axis position
    obj.mag_axis_r = equilibrium_ids_time_slice.global_quantities.magnetic_axis.r.value
    obj.mag_axis_z = equilibrium_ids_time_slice.global_quantities.magnetic_axis.z.value

    # Write the R-Z grid
    obj.grid_r = equilibrium_ids_time_slice.profiles_2d[0].grid.dim1.value
    obj.grid_dim_r = len(obj.grid_r)
    obj.grid_z = equilibrium_ids_time_slice.profiles_2d[0].grid.dim2.value
    obj.grid_dim_z = len(obj.grid_z)

    # Poloidal magnetic flux
    obj.psi = equilibrium_ids_time_slice.profiles_2d[0].psi.T / (2 * np.pi)
    obj.psi_boundary = equilibrium_ids_time_slice.global_quantities.psi_boundary / (2 * np.pi)
    obj.psi_axis = equilibrium_ids_time_slice.global_quantities.psi_axis / (2 * np.pi)
    obj.psi_sign = np.sign(obj.psi_boundary - obj.psi_axis)

    # Write the plasma current
    obj.Ip = np.abs(equilibrium_ids_time_slice.global_quantities.ip).item()

    # Instead of taking the equilibirum.time_slice[:].profiles_1d.f[-1]
    # Use the vacuum_toroidal_field.r0, b0

    obj.fpol_vacuum = (vacuum_toroidal_field_ids.b0[0] * vacuum_toroidal_field_ids.r0).item()
    # Write the FPOL
    if not equilibrium_ids_time_slice.profiles_1d.f.size:
        obj.fpol = np.asarray([obj.fpol_vacuum])
        obj.fpol_flux = np.asarray([equilibrium_ids_time_slice.global_quantities.psi_axis])
    else:
        obj.fpol = equilibrium_ids_time_slice.profiles_1d.f.value
        obj.fpol_flux = equilibrium_ids_time_slice.profiles_1d.psi.value

    # Equilibrium type
    obj.type = 'div' if equilibrium_ids_time_slice.boundary.type else 'lim'

    # Additional data to obtain
    obj.a = equilibrium_ids_time_slice.boundary.minor_radius.value
    obj.Area = equilibrium_ids_time_slice.global_quantities.area.value
    if summary_ids is not None:
        if summary_ids.global_quantities.power_loss.value.size:
            obj.Psol = np.abs(summary_ids.global_quantities.power_loss.value[0]).item()

    if correct_helicty:
        correct_equilibrium_helicity(obj)

    return obj

def createEqdskFromSlice(slice: Any, header: str="") -> EQDSKIO:
    """From a slice of the equilibrium IDS, create an EQDSKIO object.
    """
    import numpy as np
    eqObj = EQDSKIO()

    profile2d = slice.profiles_2d[0]

    # Rectangular grid
    R = profile2d.grid.dim1
    Z = profile2d.grid.dim2

    Nw = R.shape[0]
    Nh = Z.shape[0]

    eqObj.NW = Nw
    eqObj.NH = Nh

    # Setting R,Z grid
    eqObj.RLEFT = R[0]
    eqObj.RDIM = R[-1] - R[0]
    eqObj.ZMID = 0.5 * (Z[0] + Z[-1])
    eqObj.ZDIM = Z[-1] - Z[0]
    psi = profile2d.psi / (2*np.pi)
    eqObj.PSIRZ = psi.T

    profile1d = slice.profiles_1d

    def densify(f: np.ndarray, N: int) -> np.ndarray:
        """Extends the array to N if shape is different.
        """
        from scipy.interpolate import interp1d
        if f.shape[0] == 0:
            return np.array([0 for i in range(N)])

        if f.shape[0] == N:
            return f

        x = np.linspace(0, 1, f.shape[0])
        i = interp1d(x, f)
        newX = np.linspace(0, 1, N)
        newF = i(newX)
        return newF

    # print(profile1d.f.shape, profile1d.pressure.shape,
    #       profile1d.dpressure_dpsi.shape)

    eqObj.FPOL = list(densify(profile1d.f,Nw))
    eqObj.PRES = list(densify(profile1d.dpressure_dpsi, Nw) / (2 * np.pi))
    eqObj.FFPRIM = list(densify(profile1d.f_df_dpsi, Nw) / (2 * np.pi))
    eqObj.QPSI = list(densify(profile1d.q, Nw))

    globQuant = slice.global_quantities


    eqObj.RMAXIS = globQuant.magnetic_axis.r
    eqObj.ZMAXIS = globQuant.magnetic_axis.z
    eqObj.SIBRY = globQuant.psi_boundary / (2 * np.pi)
    eqObj.SIMAG = globQuant.psi_axis / (2 * np.pi)
    eqObj.CURRENT = globQuant.ip
    eqObj.RCENTR = globQuant.magnetic_axis.r
    BCENTR = profile1d.f[0] / globQuant.magnetic_axis.r
    # print(f"B_t at centre: {BCENTR} T")
    current = np.abs(globQuant.ip)
    current = current / 1e6
    current = round(current)
    # print(f"Drsep: {slice.boundary_separatrix.gap[-1].value}")
    if len(slice.boundary_separatrix.gap):
        drsep = slice.boundary_separatrix.gap[-1].value
        header = f"{current}MA {header} dsep={100 * drsep:.02f}cm - 3 {Nw} {Nh}"
    else:
        header = f"{current}MA {header} - 3 {Nw} {Nh}"
    eqObj.HEADER = header
    eqObj.BCENTR = BCENTR

    boundary = slice.boundary
    eqObj.RBBBS = boundary.outline.r
    eqObj.ZBBBS = boundary.outline.z
    eqObj.NBBBS = len(boundary.outline.r)


    return eqObj

def addWallDescriptionToEqdsk(eqObj: EQDSKIO, idsWall: Any) -> None:
    """From the wall IDS take the wall silhouette points and save it to the
    EQDSKIO object.
    """
    if len(idsWall.description_2d) == 0:
        eqObj.RLIM = [1]
        eqObj.ZLIM = [1]
        eqObj.LIMITR = 1
    else:
        description2d = idsWall.description_2d[0]
        eqObj.RLIM = description2d.limiter.unit[0].outline.r
        eqObj.ZLIM = description2d.limiter.unit[0].outline.z
        eqObj.LIMITR = len(description2d.limiter.unit[0].outline.r)

from ._iter import EquilibriumIterator

__all__ =[
    "Equilibrium",
    "EquilibriumIterator",
    "EQDSKIO",
    "EQ",
    "correct_equilibrium_helicity",
    "getEquilibriumFromIMAS",
    "getEquilibriumFromEQDSKG",
    "getEquilibriumFromIMASSlice",
    "getEquilibriumFromEQFile",
    "createEqdskFromSlice",
    "addWallDescriptionToEqdsk",
    "getOwlConlensGraph"
]