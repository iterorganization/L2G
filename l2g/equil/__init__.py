from ._equilibrium import Equilibrium, correct_equilibrium_helicity
from ._eqdskg import EQDSKIO
from ._eq import EQ

def getEquilibriumFromEQFile(eqdsk_file: str, correct_helicity: bool = True) -> Equilibrium:
    """Gather the data required from an EQDSKIO object and populate the
    Equilibrium class.
    """

    eqdsk = EQDSKIO(eqdsk_file)
    return getEquilibriumFromEQDSKG(eqdsk)

def getEquilibriumFromEQDSKG(eqdsk_obj: EQDSKIO, correct_helicity=True) -> Equilibrium:
    """Gather the data required from an EQDSKIO object and populate the
    Equilibrium class.
    """
    import numpy as np
    obj = Equilibrium()

    # Write the wall Limiter
    obj.wall_contour_r = eqdsk_obj.getRLIM()
    obj.wall_contour_z = eqdsk_obj.getZLIM()

    # Write the magnetic axis position
    obj.mag_axis_r = eqdsk_obj.getRMAXIS()
    obj.mag_axis_z = eqdsk_obj.getZMAXIS()


    # Write the R-Z grid
    RLEFT, ZMID = eqdsk_obj.getRLEFT(), eqdsk_obj.getZMID()
    RDIM, ZDIM = eqdsk_obj.getRDIM(), eqdsk_obj.getZDIM()
    NW, NH = eqdsk_obj.getNW(), eqdsk_obj.getNH()
    ZMIN = ZMID - 0.5 * ZDIM
    ZMAX = ZMID + 0.5 * ZDIM

    obj.grid_r = np.linspace(RLEFT, RLEFT + RDIM, NW)
    obj.grid_dim_r = NW

    obj.grid_z = np.linspace(ZMIN, ZMAX, NH)
    obj.grid_dim_z = NH

    # Poloidal magnetic flux
    obj.psi = np.asarray(eqdsk_obj.getPSIRZ())
    obj.psi_boundary = eqdsk_obj.getSIBRY()
    obj.psi_axis = eqdsk_obj.getSIMAG()
    obj.Ip = eqdsk_obj.getCURRENT()

    obj.psi_sign = np.sign(obj.psi_boundary - obj.psi_axis )

    # Write the Plasma current
    # Write the FPOL
    obj.fpol = np.asarray(eqdsk_obj.getFPOL())
    obj.fpol_flux = np.linspace(eqdsk_obj.getSIMAG(), eqdsk_obj.getSIBRY(), NW)
    # Write the FPOL in vacuum. Just take the last value
    obj.fpol_vacuum = obj.fpol[-1]

    # Equilibrium type
    obj.type = ''

    if correct_helicity:
        correct_equilibrium_helicity(obj)

    return obj

def getEquilibriumFromIMASSlice(shot: int, run: int, user: str = "public",
        machine: str = "iter", time: float = 0.0, data_version: str = "3",
        correct_helicity: bool = True) -> Equilibrium:

    import imas
    import imas.imasdef
    db_entry = imas.DBEntry(shot=shot, run=run, user_name=user,
        db_name=machine, backend_id=imas.imasdef.MDSPLUS_BACKEND,
        data_version=data_version)
    db_entry.open()
    summary_ids = db_entry.get("summary")
    wall_ids = db_entry.get("wall")

    equilibrium_ids = db_entry.get_slice("equilibrium", time,
                                         imas.imasdef.CLOSEST_INTERP)

    equilibrium =  getEquilibriumFromIMAS(equilibrium_ids.time_slice[0],
        wall_ids, summary_ids, correct_helicty=correct_helicity)
    return equilibrium

def getEquilibriumFromIMAS(equilibrium_ids_time_slice, wall_ids,
                           summary_ids=None, correct_helicty=True) -> Equilibrium:
    """Populate the Equilibrium class with relevant data.

    Optional:
        summary_ids can be optionally provided to obtain other relevant values.
    """
    import numpy as np
    obj = Equilibrium()

    # Write the wall limiter
    obj.wall_contour_r = wall_ids.description_2d[0].limiter.unit[0].outline.r
    obj.wall_contour_z = wall_ids.description_2d[0].limiter.unit[0].outline.z

    # Write the magnetic axis position
    obj.mag_axis_r = equilibrium_ids_time_slice.global_quantities.magnetic_axis.r
    obj.mag_axis_z = equilibrium_ids_time_slice.global_quantities.magnetic_axis.z

    # Write the R-Z grid
    obj.grid_r = equilibrium_ids_time_slice.profiles_2d[0].grid.dim1
    obj.grid_dim_r = len(obj.grid_r)
    obj.grid_z = equilibrium_ids_time_slice.profiles_2d[0].grid.dim2
    obj.grid_dim_z = len(obj.grid_z)

    # Poloidal magnetic flux
    obj.psi = equilibrium_ids_time_slice.profiles_2d[0].psi.T / (2 * np.pi)
    obj.psi_boundary = equilibrium_ids_time_slice.global_quantities.psi_boundary / (2 * np.pi)
    obj.psi_axis = equilibrium_ids_time_slice.global_quantities.psi_axis / (2 * np.pi)
    obj.psi_sign = np.sign(obj.psi_boundary - obj.psi_axis)

    # Write the plasma current
    obj.Ip = np.abs(equilibrium_ids_time_slice.global_quantities.ip)

    # Write the FPOL
    obj.fpol = equilibrium_ids_time_slice.profiles_1d.f
    obj.fpol_flux = equilibrium_ids_time_slice.profiles_1d.psi

    # Write the FPOL in vacuum. Just take the last value
    obj.fpol_vacuum = obj.fpol[-1]

    # Equilibrium type
    obj.type = 'div' if equilibrium_ids_time_slice.boundary.type else 'lim'

    # Additional data to obtain
    obj.a = equilibrium_ids_time_slice.boundary.minor_radius
    obj.Area = equilibrium_ids_time_slice.global_quantities.area
    if summary_ids is not None:
        obj.Psol = np.abs(summary_ids.global_quantities.power_loss.value[0])

    if correct_helicty:
        correct_equilibrium_helicity(obj)

    return obj

def createEqdskFromSlice(slice, HEADER="") -> EQDSKIO:
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

    eqObj.setNW(Nw)
    eqObj.setNH(Nh)

    # Setting R,Z grid
    eqObj.setRLEFT(R[0])
    eqObj.setRDIM(R[-1] - R[0])
    eqObj.setZMID((Z[0] + Z[-1]) / 2)
    eqObj.setZDIM(Z[-1] - Z[0])
    psi = profile2d.psi / (2*np.pi)
    eqObj.setPSIRZ(psi.T)

    profile1d = slice.profiles_1d

    def densify(f, N):
        """Extends the array to N if shape is different.
        """
        from scipy.interpolate import interp1d
        if f.shape[0] == 0:
            return [0 for i in range(N)]

        if f.shape[0] == N:
            return f

        x = np.linspace(0, 1, f.shape[0])
        i = interp1d(x, f)
        newX = np.linspace(0, 1, N)
        newF = i(newX)
        return newF

    # print(profile1d.f.shape, profile1d.pressure.shape,
    #       profile1d.dpressure_dpsi.shape)
    FPOL = densify(profile1d.f,Nw)

    eqObj.setFPOL(FPOL)
    eqObj.setPRES(densify(profile1d.pressure, Nw))
    eqObj.setPPRIME(densify(profile1d.dpressure_dpsi, Nw) / (2 * np.pi))
    eqObj.setFFPRIM(densify(profile1d.f_df_dpsi, Nw) / (2 * np.pi))
    eqObj.setQPSI(densify(profile1d.q, Nw)) # empty

    globQuant = slice.global_quantities


    eqObj.setRMAXIS(globQuant.magnetic_axis.r)
    eqObj.setZMAXIS(globQuant.magnetic_axis.z)
    eqObj.setSIBRY(globQuant.psi_boundary / (2 * np.pi))
    eqObj.setSIMAG(globQuant.psi_axis / (2 * np.pi))
    eqObj.setCURRENT(globQuant.ip)
    eqObj.setRCENTR(globQuant.magnetic_axis.r)
    BCENTR = profile1d.f[0] / globQuant.magnetic_axis.r
    # print(f"B_t at centre: {BCENTR} T")
    current = np.abs(globQuant.ip)
    current = current / 1e6
    current = round(current)
    # print(f"Drsep: {slice.boundary_separatrix.gap[-1].value}")
    if len(slice.boundary_separatrix.gap):
        drsep = slice.boundary_separatrix.gap[-1].value
        HEADER = f"{current}MA {HEADER} dsep={100 * drsep:.02f}cm - 3 {Nw} {Nh}"
    else:
        HEADER = f"{current}MA {HEADER} - 3 {Nw} {Nh}"
    eqObj.setHEADER(HEADER)
    eqObj.setBCENTR(BCENTR)

    boundary = slice.boundary
    eqObj.setRBBBS(boundary.outline.r)
    eqObj.setZBBBS(boundary.outline.z)
    eqObj.setNBBBS(len(boundary.outline.r))


    return eqObj

def addWallDescriptionToEqdsk(eqObj, idsWall):
    """From the wall IDS take the wall silhouette points and save it to the
    EQDSKIO object.
    """
    if len(idsWall.description_2d) == 0:
        eqObj.setRLIM([1])
        eqObj.setZLIM([1])
        eqObj.setLIMITR(1)
    else:
        description2d = idsWall.description_2d[0]
        eqObj.setRLIM(description2d.limiter.unit[0].outline.r)
        eqObj.setZLIM(description2d.limiter.unit[0].outline.z)
        eqObj.setLIMITR(len(description2d.limiter.unit[0].outline.r))

from ._iter import EquilibriumIterator