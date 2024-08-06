r"""Contains methods for the ELM PLM method. Look at 'A model of ELM filament
energy evolution due to parallel losses' by Fundamenski
(doi:10.1088/0741-3335/48/1/008), section 5.

One of the main assumptions is the usage of the "Parallel dynamics with
acoustic response", section 5.2, where:

    M = 1; Mach number
    \alpha_i = 3.5 \xi; ion flux limiter
    \alpha_e = 5 \xi (v_ti / v_te); electron flux limiter

The method solves the equations 5.5 in dimensionless form.

"""

from scipy.interpolate import interp1d
import logging
import numpy as np
import math
log = logging.getLogger(__name__)


def tau_ee() -> float:
    r"""Calculate the factor of the Tau_{e-e} the electron-electron collision
    time.

    \tau_{e-e} = \frac{12 \pi^{3/2} \eps_0^2 m_e^{1/2}}{2^{1/2} e^4 ln(\Lambda)} \frac{T_e^{3/2}}{n_e}
    """
    pi = np.pi
    eps0 = 8.854187 * 1e-12 #  F/m
    me = 9.109e-31 # kg
    ln_delta = 20 # Couloumb constant
    e = 1.602176634e-19 # C


    c = (12 * pi**(3/2) * eps0**2 * np.sqrt(me)) / (np.sqrt(2) * e**4 * ln_delta)

    c *= e**(3/2) # This comes from the T_e^{3/2} which is in electron volts.
    return c

def tau_ii() -> float:
    r"""Calculate the factor of the Tau_{i-i} the ion-ion collision time.

    \tau_{i-i} = \frac{12 \pi^{3/2} \eps_0^2 m_p^{1/2}}{Z^4 2^{1/2} e^4 ln(\Lambda)} \frac{T_e^{3/2}}{n_e}
    """
    pi = np.pi
    eps0 = 8.854187e-12 #  F/m
    A = 2 # Atomic mass average for D-D plasmas. 2.5 for D-T
    mp = 1.67262192e-27 # kg
    ln_delta = 20 # Couloumb constant
    e = 1.602176634e-19 # C
    A = 1

    c = (12 * pi**(3/2) * eps0**2 * np.sqrt(A*mp)) / (np.sqrt(2) * e**4 * ln_delta)

    c *= e**(3/2)
    return c

def tau_ei() -> float:
    r"""Calculate the factor of the Tau_{e-i} the electron-ion collision
    time.

    \tau_{e-e} = \frac{12 \pi^{3/2} \eps_0^2 m_e^{1/2}}{2^{1/2} e^4 ln(\Lambda)} \frac{T_e^{3/2}}{n_e}
    """
    pi = np.pi
    eps0 = 8.854187 * 1e-12 #  F/m
    me = 9.109e-31 # kg
    mp = 1.67262192e-27 # kg

    ln_delta = 20 # Couloumb constant
    e = 1.602176634e-19 # C

    A = 1
    c = (12 * pi**(3/2) * eps0**2 * np.sqrt(me)) / (np.sqrt(2) * e**4 * ln_delta)

    c *= e**(3/2) # This comes from the T_e^{3/2} which is in electron volts.
    return c

def tau_ie() -> float:
    r"""Calculate the factor of the Tau_{i-e} the ion-electron collision
    time.

    \tau_{e-e} = \frac{12 \pi^{3/2} \eps_0^2 m_i}{2^{1/2} m_e^{1/2} e^4 ln(\Lambda)} \frac{T_e^{3/2}}{n_e}
    """
    pi = np.pi
    eps0 = 8.854187 * 1e-12 #  F/m
    me = 9.109e-31 # kg
    mp = 1.67262192e-27 # kg

    ln_delta = 20 # Couloumb constant
    e = 1.602176634e-19 # C

    A = 1
    c = (6 * pi**(3/2) * eps0**2 * A * mp) / (np.sqrt(2) * np.sqrt(me)  * e**4 * ln_delta)

    c *= e**(3/2) # This comes from the T_e^{3/2} which is in electron volts.
    return c

class ELM_PLM(object):
    def __init__(self, drsep, conlen):
        self.drsep = drsep
        self.conlen = conlen
        self.interp_conlen = interp1d(drsep, conlen, kind="linear",
                                      fill_value="extrapolate")

        ## Calculated and populated with calculate_loss_profiles
        self.Ei0: float = None
        self.ei: np.ndarray = None # Normalized ion energy

        self.Ee0: float = None
        self.ee: np.ndarray= None # Normalized electron energy

        self.t: np.ndarray = None # Parametric time

        self.tau_n0: float = None # Number density loss time. Also used to normalize
                           # parametric time
        self.n: np.ndarray = None # Normalized number density
        self.n0: float = None

    def calculate_loss_profiles(self, v_elm=500.0, n0 = 0.375, Te0=2500, Ti0=2500,
            Aion=2.5, Z=1, r_max=0.5):
        """Calculates the profiles of the average density and energy of the
        ELM filament.

        Fundamenski (doi:10.1088/0741-3335/48/1/008), section 5.

        Arguments:
            v_elm (float): Radial propagation speed [m/s].
            n0 (float): Filament peak density at R=Rsep [10^20 m^-3].
            Te0 (float): Filament peak electron temperature at R=Rsep [eV]
            Ti0 (float): Filament peak ion temperature at R=Rsep [eV]
            Aion (float): Average atomic mass. 2.5 for DT plasma.
        """

        MP = 1.673E-27 # Mass of a proton (kg)
        ME = 9.108E-31 # Mass of an electron (kg)
        QE = 1.602E-19 # Charge of an electron (Coulombs)

        # MACH=1.0 based on 5.2 Parallel dynamics with acoustic responses
        MACH = 1.0

        # Arrays
        a_n = []
        a_Te = []
        a_Ti = []
        a_ee = []
        a_ei = []
        a_W = []
        a_t = []
        a_M = []

        ## Calculate the initial density loss time in SI units
        xi0 = (1 + Z / (Ti0/Te0))**0.5
        v_ti0 = ((Ti0*QE)/(MP*Aion))**0.5 # Initial ion thermal velocity
        # Warm ion sound speed c_s = xi0 * v_ti9
        # Density loss time
        tau_n0 = self.interp_conlen(0.0) / (MACH * xi0 * v_ti0)

        ## Initial average energy of the filamen in MJ
        Ee0 = (3/2) * n0 * 1e20 * Te0 * QE
        Ei0 = (3/2) * n0 * 1e20 * Ti0 * QE

        # CONSTANTS = [n0, Ee0, Ei0, tau_n0, Aion, MP, ME, Z, QE, self.interp_conlen]
        def deriv(t:float, y: np.ndarray)->np.ndarray:
            """Derivative function that modified the yp array!
            """
            # Current dimensionless values at (r, t)
            n, ee, ei = y
            r = t * tau_n0 * v_elm
            # Temperatures in eV
            Te = (2.0/3.0) * (ee/n) * (Ee0/(n0*1e20*QE))
            Ti = (2.0/3.0) * (ei/n) * (Ei0/(n0*1e20*QE))
            # Thermal speeds
            vti = (Ti*QE / (MP*Aion))**0.5
            vte = (Te*QE / ME)**0.5
            # Connection length at this position r
            L = self.interp_conlen(r) # Get the connection length
            MACH = 1 # Constant
            Xi = (1 + Z / (Ti / Te))**0.5 # Ion speed of sound factor Xi

            ## Particle collision times.
            # E-e collision time in s.
            Tau_ee = 1.720353e-10 * (Te**(1.5)) / (n*n0)
            # I-i collision time in s. Same thing here. The temperature is in
            # eV and the density is in units of 1e20 m^-3...
            Tau_ii = 7.371932e-09 * (Aion**0.5)*(Z**-3)*(Ti**(1.5)) / (n*n0)
            # I-e (Ion-electron!) collision time in s.
            Tau_ie = 1.579482e-07 * Aion*(Z**-2)*(Te**(1.5)) / (n*n0)

            # Sonic speeds
            c_par_e = 3.2
            c_par_i = 3.9
            ## Flux limiters
            alpha_i = 3.5 * Xi
            alpha_e = 5 * Xi * (vti/vte)
            ## Calculate the time losses connected to particles and energies
            # Particle loss time (s)
            Tau_n = L / (MACH * Xi * vti)

            # Now calculate the NORMALIZED energy time losses.
            Tau_Chi_i = 2.5 * L**2 * (1 + (c_par_i*vti*Tau_ii)/(alpha_i*L)) / (c_par_i * vti**2 * Tau_ii)
            Tau_Chi_e = 2.5 * L**2 * (1 + (c_par_e*vte*Tau_ee)/(alpha_e*L)) / (c_par_e * vte**2 * Tau_ee)

            Tau_E_e = 0.6 * (1/Tau_n  + 1/Tau_Chi_e)**-1
            Tau_E_i = 0.6 * (1/Tau_n  + 1/Tau_Chi_i)**-1
            yp = np.empty(3)
            yp[0] = -y[0] / Tau_n
            yp[1] = -y[1] / Tau_E_e + (y[2] - y[1]) / Tau_ie
            yp[2] = -y[2] / Tau_E_i - (y[2] - y[1]) / Tau_ie
            yp *= tau_n0 # Normalize the equation to the dimensionless form
            return yp
        ## Initial time step

        tstep = (np.min(np.diff(self.drsep)) / v_elm)/ tau_n0
        # tstep = 0.3 * 1e-4 / (tau_n0*v_elm)
        # tstep = 1e-6
        ## Initial Condition vector y (y = [n_, Ee_, Ei_]), dimensionless
        t = 0
        y = np.array([1, 1, 1])
        with np.errstate(invalid="raise"):
            yp = deriv(0, y)
        a_t.append(t)
        a_n.append(y[0])
        a_ee.append(y[1])
        a_ei.append(y[2])

        k2 = np.empty(3)
        k3 = np.empty(3)
        k4 = np.empty(3)
        k5 = np.empty(3)
        k6 = np.empty(3)
        sol = np.empty(3)

        h = tstep

        r = 0.0

        # It's important to check when the r get's almost close to r_max
        while (not np.allclose(r, r_max)) and t < 5:
            # Do the RKF 45 calculation
            while 1:
                # Progress the computation
                # Calculate at k2
                with np.errstate(invalid="raise"):
                    try:
                        ch = h / 4
                        k6 = y + ch * yp
                        k2 = deriv(t+h, k6)

                        # Calculate k3
                        ch = 3.0 * h / 32.0
                        k6 = y + ch * (yp + 3.0 * k2)

                        # Calculate the k3 parameter
                        k3 = deriv(t + 3.0 * h / 8.0, k6)

                        # Calculate k4
                        ch = h / 2197.0;
                        k6 = y + ch * (1932.0 * yp  + (7296.0 * k3 - 7200.0 * k2))
                        k4 = deriv(t + 12.0 * h / 13.0, k6)

                        # Calculate k5
                        ch = h / 4104.0
                        k6 = y + ch * ((8341.0 * yp - 845.0 * k4) + (29440.0 * k3 - 32832.0 * k2))
                        k5 = deriv(t + h, k6)

                        # Calculate k6
                        ch = h / 20520.0
                        k2 = y + ch * ((-6080.0 * yp + (9295.0 * k4 - 5643.0 * k5)) + (41040 * k2 - 28352.0 * k3))
                        k6 = deriv(t + h / 2.0, k2)

                        # Calculate the solution at t + h
                        ch = h / 7618050.0
                        sol = y + ch * ((902880.0 * yp + (3855735.0 * k4 - 1371249.0 * k5)) + (3953664.0 * k3 + 277020.0 * k6))

                        et = np.abs(y) + np.abs(sol) + 2
                        ee = np.abs((-2090.0 * yp + (21970.0 * k4 - 15048.0 * k5)) + (22528.0 * k3 - 27360.0 * k6))
                        esttol = np.abs(h) * np.max(ee / et) / 752400.0

                        # Try to also calculate the derivative on solution to
                        # trigger FloatingPointError
                        deriv(t+h, sol)
                    except FloatingPointError:
                        # Exception raised due to division with 0.
                        # Halving the scale.
                        h *= 0.1
                        continue
                if esttol <= 1.0:
                    # Successful step!
                    break
                if esttol < 59049.0:
                    scale = 0.9 / esttol**(0.2)
                else:
                    scale = 0.1

                if h * scale < 1e-7:
                    h = 1e-7
                else:
                    h *= scale
            # Update the new values
            y = np.copy(sol)
            t += h
            yp = deriv(t, y)
            # The problem is super unstable so do not activate this to make
            # it "faster"!
            # if esttol > 0.0001889568:
            #     scale = 0.9 / esttol**0.2
            # else:
            #     scale = 5.0
            # if step_failed:
            #     scale = min(scale, 1.0)
            # h *= scale

            a_t.append(t)
            a_n.append(y[0])
            a_ee.append(y[1])
            a_ei.append(y[2])
            r = t * v_elm * tau_n0

        a_t = np.asarray(a_t)
        a_n = np.asarray(a_n)
        a_ee = np.asarray(a_ee)
        a_ei = np.asarray(a_ei)

        # Store results
        self.Ei0 = Ei0
        self.Ti = Ti0 * a_ei
        self.ei = a_ei

        self.Ee0 = Ee0
        self.Te = Te0 * a_ee
        self.ee = a_ee


        self.tau_n0 = tau_n0
        self.t = a_t
        self.n = a_n
        self.n0 = n0
        self.r = a_t * v_elm * tau_n0

    def calculate_heat_load_profile(self, E_total: float=0.6e6, nfil: int=10,
            sigmar: float=3e-2, sigmaz: float=0.45, v_elm:float = 500.0,
            fELM:float=33, rg2f: float=2):
        """Calculates the heat load of the ELM loss profiles.

        https://doi.org/10.1088/1741-4326/acdf02. Section 3, equations (2),
        (3), (4).

        Arguments:
            E_total (float): Total energy released per ELM
            nfil (int): Number of filaments.
            sigmar (float): Filament radial width [m]
            sigmaz (float): Filament vertical width [m]
            v_elm (float): Radial velocity of the ELM.
            fELM (float): ELM frequency.
            rg2f (float): Gap-to-filament ratio

        """

        if self.n is None:
            return

        t = self.t * self.tau_n0
        r = self.r

        ## 1) Energy of the filament!
        # Calculate the normalized dimensionless energy of the filament
        E = 0.5*(self.ee + self.ei)
        ## 2) Get the q|| averaged through time profile.

        # How the energy of the filament drops over time.
        dWdt = np.zeros(self.t.shape)
        dWdt[0] = (E[1]-E[0])/(t[1]-t[0]) # Forward finite difference
        dWdt[1:-1] = (E[2:] - E[:-2]) / (t[2:] - t[:-2]) # Central finite difference
        dWdt[-1] = (E[-1]-E[-2])/(t[-1] - t[-2]) # Backward finite difference

        # Get the parallel energy of a filament averaged over time
        Epar = np.zeros(r.shape)

        # Calculate the Qpar, Equation 4 from doi.org/10.1088/1741-4326/acdf02

        # First the energy per filament on every radial point. Which is
        # integrated over time to get the average value.
        C = -0.5 * dWdt / (2*np.pi*sigmar*sigmaz)
        diff_t = np.diff(t)
        for i in range(r.shape[0]):
            # Equation 2
            # r here (not r[i]) is the same as tv_r since r is derived as
            # r = t tau_n0 v_elm
            q_par_fil = C * np.exp(-(r[i] - r)**2/(2*sigmar**2))
            # But integration takes over time!
            Epar[i] = np.sum(0.5*(q_par_fil[1:] + q_par_fil[:-1])*diff_t)

        # Take into account that the parallel energy of the elm is normalized,
        # therefore in order to get the total contribution of all ELMs, we
        # distribute the total energy, divide it by the number of filaments and
        # then multiply by the frequency.

        # NB: This needs to be re-evaluated. In other words if the dWdT is
        # showing how a SINGLE filament is losing power, then in order to get
        # the total ELM contribution to be used as the input profile we need
        # to actually integrate all of the filaments that occur in a single
        # second. Otherwise use the following expression to get the ELM
        # contribution.

        # rg2f - gap to filament ratio already involves the pitch, so the
        # output is trully the q||, perpendicular to the midplane. No need to
        # further multiply it to obtain the correctly pitched profile on the
        # midplane. For consistency with the previous implementation this
        # term is still used, instead of the actual definition of the term.
        qpar = Epar * (E_total / nfil) * fELM / (rg2f)


        # Characterstic e-folding length of WELM
        efold = -((E[1:]-E[:-1])/(r[1:]-r[:-1]) / (0.5*(E[1:]+E[:-1])))**-1
        # Save the values
        self.qpar = qpar
        self.efold = efold

    def create_elm_plm_graphs(self, xlim=(-2.5, 400.0),
            save_to_file: bool=False, file_path: str=""):
        import matplotlib.pyplot as plt
        f, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True, dpi=100)

        r = self.r * 1000
        XLABEL = "$r-r_{sep}$ [mm]"
        axs[0].set_title('a) $L_{conn}[m]$')
        axs[1].set_title('b) Normalized n, $\epsilon_e$, $\epsilon_i$')

        for ax in axs:
            ax.grid(True)
            ax.set_xlabel(XLABEL)
            ax.set_xlim(xlim)

        axs[0].semilogy(self.drsep * 1e3, self.conlen)
        r = self.r * 1e3
        axs[1].semilogy(r, self.n, label="n")
        axs[1].semilogy(r, self.ee, label=r"$\epsilon_e$")
        axs[1].semilogy(r, self.ei, label=r"$\epsilon_i$")
        axs[1].legend()

        if save_to_file:
            f.savefig(file_path)
            f.clf()
            plt.close(f)
            del f
            return
        return f, axs

    def create_graphs(self, xlim=(-2.5, 400.0), save_to_file: bool=False,
            file_path: str=""):
        import matplotlib.pyplot as plt
        f, axs = plt.subplots(1, 5, figsize=(12, 6), tight_layout=True, dpi=100)

        r = self.r * 1000
        XLABEL = "$r-r_{sep}$ [mm]"
        axs[0].set_title('a) $L_{conn}[m]$')
        axs[1].set_title('b) $W_{fil}/W_{0fil}$')
        axs[2].set_title(r'c) $\lambda [mm]$')  # $[(dW/dr)/W]^{-1}[mm]$
        axs[3].set_title(r'd) $q_{\parallel, ELM}$ $[\frac{MW}{m^2}$]')
        axs[4].set_title(r'e) $T_e, T_i$ [eV]')


        for ax in axs:
            ax.grid(True)
            ax.set_xlabel(XLABEL)
            ax.set_xlim(xlim)

        axs[0].semilogy(self.drsep * 1e3, self.conlen)

        axs[1].semilogy(r,0.5*(self.ee + self.ei), 'b',
                      label='Steady state')

        axs[2].semilogy(r[:-1], self.efold * 1e3, 'b',)

        axs[3].semilogy(r, self.qpar*1e-6, 'y--')

        axs[4].semilogy(r, self.Te, 'y--', label="Te")
        axs[4].semilogy(r, self.Ti, 'r--', label="Ti")
        axs[4].legend()

        if save_to_file:
            f.savefig(file_path)
            f.clf()
            plt.close(f)
            del f
            return
        return f, axs

def test():


    # For the Parallel Loss model!
    Aion = 2.5  # 2.5 for DT plasma
    n0 = 0.375  # filament peak density at R=Rsep [10^20 m^-3]
    Te0 = 2500  # filament peak electron temperature at R=Rsep [eV]
    Ti0 = Te0  # filament peak ion temperature at R=Rsep [eV]
    v_elm = 500  # radial propagation speed [m/s]
    # For the calculation of heat load profile
    fELM =33
    nfil=6
    poloWidth=0.6
    aveFac=0.1


    import numpy as np
    import l2g.equil
    import l2g.hlm.steady_state
    import matplotlib.pyplot as plt


    eq_file = "shot_105045_run_1_t_75.028s.eqdsk"
    eq_file = "../equil/g900003_00230_ITER_15MA_eqdsk16HR.txt"

    equilibrium = l2g.equil.getEquilibriumFromEQFile(eq_file)

    # Equilibrium analysis object
    eq = l2g.external.equilibrium_analysis.EQA()
    eq.setEquilibrium(equilibrium)
    eq.evaluate()
    print("getting owl connection length graph")
    drsep, conlen_down, conlen_up = l2g.equil.getOwlConlensGraph(eq)
    print(drsep.size, conlen_down.size)

    # See where the connection lengths drop to zero.
    zero_drop_conlen_down = np.where(np.isclose(conlen_down, 0))[0]
    if len(zero_drop_conlen_down):
        r_last_index = zero_drop_conlen_down[0]
    else:
        r_last_index = -1

    drsep = drsep[:r_last_index]
    conlen_down = conlen_down[:r_last_index]
    conlen_up = conlen_up[:r_last_index]


    r_out = drsep
    # eq_file = "../../eqdsks/shot_105045_run_1_t_75.028s.eqdsk"
    # equilibrium = l2g.equil.getEquilibriumFromEQFile(eq_file)
    # eq = l2g.equil.EQ(equilibrium)
    # eq.evaluate()



    elm_loss = l2g.hlm.steady_state.ELM_loss(v_elm, n0, Te0, Ti0, Aion,
                            fELM, nfil, poloWidth, aveFac, conlen_down, r_out)

    [n, Te, Ti, Ee, Ei,
     W, t_out, M] = elm_loss.ELM_loss_profile(v_elm, n0, Te0, Ti0, r_out,
                                              conlen_down)

    f = plt.figure()
    ax = f.add_subplot()
    # ax.plot(t_out.NORM, Ei.NORM/Ee.NORM)
    ax.plot(t_out.NORM, Ei.NORM/Ee.NORM)

    f = plt.figure()
    ax = f.add_subplot()
    # ax.plot(t_out.NORM, n.NORM, label="n")
    # ax.plot(t_out.NORM, Ee.NORM, label="Ee")
    # ax.plot(t_out.NORM, Ei.NORM, label="Ei")
    ax.set_yscale("log")
    ax.set_title("orig")
    ax.set_ylim((1e-7,1))
    # ax.set_xlim((0.0, 0.003))
    ax.set_xlabel(r"$\Delta_{sep}$ - radial distance along the midplane [m$")
    ax.legend()


    obj = ELM_PLM(drsep, conlen_down)
    obj.calculate_loss_profiles(v_elm=v_elm, n0=n0, Te0=Te0, Ti0=Ti0,
                                Aion=Aion, Z=1, r_max=drsep[-1])

    a_t = obj.t
    a_r = obj.r
    a_n = obj.n
    a_Ee = obj.ee
    a_Ei = obj.ei

    ax.plot(a_r, a_n, label="n")
    ax.plot(a_r, a_Ee, label="Ee")
    ax.plot(a_r, a_Ei, label="Ei")
    ax.set_ylim((1e-7,1))
    # ax.set_xlim((0.0, 0.003))
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta_{sep}$ - radial distance along the midplane [m]")
    ax.set_title("Results of PLM profiles for flat-top AFP (s/r: 105045/1, t:75.028s)")
    ax.legend()
    # plt.show()

    f = plt.figure()
    ax = f.add_subplot()
    a_Ei = np.array(a_Ei)
    a_Ee = np.array(a_Ee)
    ax.plot(a_t, a_Ei/a_Ee)
    ax.set_title("new")

    # f = plt.figure()
    # ax = f.add_subplot()
    # ax.plot(data[:, 0], data[:, 3])
    # ax.set_yscale("log")
    plt.show()



    # Now for the powers
    WELM0 = 0.6  # Energy lost per ELM [MJ]
    fELM = 33 # ELM frequency [Hz]
    nfil = 10  # number of filaments
    poloWidth = 0.6  # poloidal width of filament [m]
    Etot = 0.6e6  # prescribed total energy released per ELM, dWELM [W]
    aveFac = 1.7  # averaging factor
    sigmar = 3e-2  # filament radial width [m]
    sigmaz = 0.45  # filament vertical height [m]
    lambdaW = 0.0213  # approximate radial e-folding length of the ELM
                      # energy (fit for drsep1 = 0-10 cm for the QDN
                      # equilibrium)


    data_jnme = elm_loss.calculate_ELM_profile(W, r_out, WELM0, fELM, nfil,
        poloWidth, aveFac, t_out, v_elm, Etot, sigmar, sigmaz, lambdaW,
        npts=len(r_out), rlim=r_out[-1])

    obj.calculate_heat_load_profile(nfil=10)

    ff, axs = plt.subplots(1, 4, figsize=(12, 6), tight_layout=True, dpi=100)
    axs[0].set_title('a) $L_{conn}[m]$')
    axs[1].set_title('b) $W_{fil}/W_{0fil}$')
    axs[2].set_title('c) $\lambda [mm]$')  # $[(dW/dr)/W]^{-1}[mm]$
    axs[3].set_title(r'a) $q_{\parallel, ELM} + q_{\parallel, inter-ELM}[\frac{MW}{m^2}]$')
    XLABEL = '$r-r_{sep}[mm]$'
    axs[0].semilogy(r_out * 1e3, conlen_down)
    axs[0].grid(True)
    axs[0].set_xlabel(XLABEL)
    axs[0].set_xlim((-2.5, 200.0))

    axs[1].grid(True)
    axs[1].set_xlim((0, 200.0))
    axs[1].set_xlabel(XLABEL)
    axs[1].semilogy(r_out[:-1]*1000, W.NORM[:-1], 'b',
                  label='Steady state')

    axs[2].semilogy(np.array(data_jnme[0][:-1]) * 1000,
                    np.array(data_jnme[2]) * 1000, 'b',label="OLD")
    axs[2].semilogy(v_elm * obj.tau_n0*obj.t[:-1] * 1000,
                    obj.efold * 1000, 'r-',label="NEW")
    axs[2].legend()
    axs[2].grid(True)
    axs[2].set_xlabel(XLABEL)
    axs[2].set_xlim((-2.5, 200.0))
    print("data_jnme[1]", data_jnme[1])
    axs[3].semilogy(data_jnme[0] * 1000, data_jnme[1], 'y--',
                    label='ELM PLM')
    axs[3].semilogy(v_elm * obj.tau_n0 * obj.t * 1000, obj.qpar*1e-6, 'k--',
                    label='NEW')
    axs[3].legend()
    axs[3].set_xlabel(XLABEL)
    axs[3].set_xlim((-2.5, 400.0))
    plt.show()

    obj.create_elm_plm_graphs()
    plt.show()

    obj.create_graphs()
    plt.show()

    import l2g.equil
    x = l2g.equil.getEquilibriumFromEQFile("shot_105045_run_1_t_75.028s.eqdsk")
    eq = l2g.equil.EQ(x)
    eq.evaluate()
    Rb, Z, Btotal, Bpm = eq.get_midplane_info(which="owl")

    x = v_elm * obj.tau_n0 * obj.t
    area = np.pi*((Rb+x[-1])**2 - Rb**2)
    print(area)
    # area = 1
    qpar = obj.qpar * 1e-6 #* Btotal/Bpm * sigmaz/(np.pi *Bpm * Rb)
    print("gap2fil", np.pi * Rb * Bpm / (sigmaz * Btotal * nfil))
    qinterp_new = interp1d(x, qpar)


    def fun_new(x):
        return 2 * 2 * np.pi * (Rb + x) * qinterp_new(x)

    def fun_old(x):
        return 2 * 2 * np.pi * (Rb + x) * qinterp_old(x)

    from scipy.integrate import quad
    I = quad(fun_new, 0, x[-1])[0]
    E = Etot * fELM * 1e-6
    print(I, E, I/E)


if __name__ == "__main__":
    test()