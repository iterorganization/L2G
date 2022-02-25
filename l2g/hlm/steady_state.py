import genericpath
import numpy as np
import math
from scipy.interpolate import interp1d
import logging
log = logging.getLogger(__name__)

def conservative(drsep: np.ndarray) -> np.ndarray:
    """From input drsep (in this case it has to be distance from the second
    separatrix) calculate the contributions from inter-ELM and ELMs (old HNLS).

    Heat load parallels are hard-coded, so are decay lengths!

    .. math::

       q_{elm \parallel} & = 5 \frac{MW}{m^2} \\\\
       \lambda_{q, elm} & = 0.09 m \\\\
       q_{inter-elm \parallel} & = 3 \frac{MW}{m^2} \\\\
       \lambda_{q, inter-elm} & = 0.17 m
    """
    q_elm = 5e6 * np.exp(-(drsep) / 0.09)
    q_inter = 3e6 * np.exp(-(drsep) / 0.17)
    return q_elm + q_inter

def inter_ELM(drsep: np.ndarray, R_bdry: float, B_total: float,
              B_pol: float, P_sol: float = 100e6, Rb: float = 0.025,
              lambda_n: float = 0.005, lambda_m: float = 0.17,
              F: float = 0.5) -> np.ndarray:
    """Calculates the q_parallel heat load profile using the update inter-ELM
    plasma profile for inter-ELM contribution (Steady-State).

    Returns:
        q (arr): Inter ELM parallel heat load profile.

    Arguments:
        drsep (arr): 1D array of distances from the midplane. In meters
        R_bdry (float): Starting radial boundary of midplane.
        B_total (float): Total magnetic magnitude at midplane start.
        B_pol (float): Poloidal component of magnetic field at midplane start.
        Rb (float): Breakdown. Default 25 mm.
        lambda_n (float): Near decay length. Default 5 mm.
        lambda_m (float): Main decay length. Default 170 mm.
        P_sol (float): Power at SOL. Default 100 MW.
        F (float): Power split between midplanes. Default 0.5.
    """
    if not(isinstance(drsep, np.ndarray)):
        if not(isinstance(drsep, list)):
            drsep = [drsep]
        drsep = np.array(drsep, dtype=np.float64)

    q = np.zeros(drsep.shape)

    C = B_total * F * P_sol / (2 * np.pi * R_bdry * B_pol * lambda_n)

    q = np.where(drsep < Rb,
                 C * np.exp(-drsep / lambda_n),
                 C * np.exp(-drsep / lambda_m + Rb / lambda_m - Rb / lambda_n))
    return q

def ELM(drsep: np.ndarray, elm_data_r: np.ndarray,
        elm_data_q: np.ndarray) -> np.ndarray:
    """Returns ELM contributed q_parallal heat load profile. This function
    itself does not evaluate the ELM profile, this is performed by a separate
    command (see get_elm_data).

    Returns:
        q (arr): ELM heat load profile

    Arguments:
        drsep (arr): 1D array of distances from the midplane. In meters.
        elm_data_r (arr): 1D array of X-axis of the numerical ELM profile.
            In meters.
        elm_data_q (arr): 1D array of numerical values of ELM. In Watts.
    """
    if not(isinstance(drsep, np.ndarray)):
        if not(isinstance(drsep, list)):
            drsep = [drsep]
        drsep = np.array(drsep, dtype=np.float64)

    interp = interp1d(elm_data_r, elm_data_q, fill_value='extrapolate')
    q = interp(drsep)
    return q

class UNN():
    '''Documentation for UNN

    '''
    def __init__(self, N=0):
        super(UNN, self).__init__()
        self.NORM = np.zeros(N)
        self.UNNORM = np.zeros(N)

    def setNORM(self, a):
        self.NORM = a

    def setUNNORM(self, a):
        self.UNNORM = a

    def getNORM(self):
        return self.NORM

    def getUNNORM(self):
        return self.UNNORM


class ELM_loss(object):
    """Documentation for ELM_loss

    """
    def __init__(self, v_elm, n0, Te0, Ti0, Aion, fELM, nfil, poloWidth,
                 aveFac, Lcon, r):
        super(ELM_loss, self).__init__()
        self.v_elm = v_elm
        self.n0 = n0
        self.Te0 = Te0
        self.Ti0 = Ti0
        self.Aion = Aion
        self.fELM = fELM
        # self.WELM = WELM
        self.nfil = nfil
        self.L = Lcon

        # Interpolate the Connection length data as a function of drsep r

        self.interp1d_conLen = interp1d(r, Lcon, fill_value='extrapolate')

    def get_M(self, t):
        '''Function to calculate the mach number as a function of time. This
        is modifiable according the the SOL regime we're in. See sections 5.1,
        5.2 and 5.3 of Wojtek's PPCF paper.
        '''
        return 1

    def RK4(self, r, tin, yin, h, fun, funArgs):
        '''Function to update the current vector of values y using the
        Runge-Kutta 4th order solver. The solver updates the system of
        equations dy/dt=f(t,y) where y is a vector and f is a vector of
        functions (called fun here). The h is the step size.
        '''
        # Start by evaluating k1=f(yin,tin), the slope at the beginning of
        # the interval
        k1 = fun(r, tin, yin, funArgs)
        # Start by evaluating k1=f(yin,tin), the slope at the beginning of
        # the interval
        k1 = fun(r, tin, yin, funArgs)
        # Now the slope at the midpoint of the interval, extrapolating using k1
        k2 = fun(r, tin+h/2, yin+h/2*k1, funArgs)
        # And now the slope at the midpoint of the interval,
        # extrapolating using k2
        k3 = fun(r, tin+h/2, yin+h/2*k2, funArgs)
        # Finally the slope at the end of the interval, extrapolating using k3
        k4 = fun(r, tin+h, yin+h*k3, funArgs)
        # RK4 method uses y_(n+1)=y_n+h/6*(k1+2*k2+2*k3+k4):
        yout = yin+h/6*(k1+2*k2+2*k3+k4)
        return yout

    def f(self, r, t, y, args):
        """Function to calculate dy/dt' according to equations (5.5), where
        y=[n',Ee',Ei'].
        """
        Ee0, Ei0, QE, n0, A, MP, ME, Z, tau_n0 = args


        n_ = y[0]
        Ee_ = y[1]
        Ei_ = y[2]
        # global Ee0
        # global Ei0
        # global QE
        # global n0
        # global A
        # global MP
        # global ME
        # global Z
        # global tau_n0
        L = self.get_L(r)  # Calculate the connection length at this time
        M = self.get_M(t)  # Calculate the mach number at this time
        n = n_*self.n0  # Density in SI units
        Te = (2/3)*((Ee_*Ee0)/(n*1E20*QE))  # Electron temperature in eV
        Ti = (2/3)*((Ei_*Ei0)/(n*1E20*QE))  # Ion temperature in eV
        xi = (1+Z/(Ti/Te))**0.5
        v_ti = ((Ti*QE)/(MP*A))**0.5  # ion thermal speed
        v_te = ((Te*QE)/ME)**0.5  # electron thermal speed
        tau_n = L/(M*xi*v_ti)  # Particle loss time (s)
        tau_ee = 1.67E-10*((Te)**(3/2))*(1/n)  # e-e collision time in s
        tau_ii = 7.1E-9*(A**0.5)*(Z**-3)*((Ti)**(3/2))*(1/n)  # i-i
                                                              # collision
                                                              # time in s
        [alpha_e, alpha_i] = self.flx_lim(M, xi, v_ti, v_te)  # Calculate the
                                                         # flux limiters
        #  Electron temperature loss time (s):
        tau_chi_e = (5/2)*(L**2)*((1+((3.2*v_te*tau_ee) /
                                      (alpha_e*L)))/(3.2*(v_te**2)*tau_ee))
        #  Ion temperature loss time (s):
        tau_chi_i = (5/2)*(L**2)*((1+((3.9*v_ti*tau_ii) /
                                      (alpha_i*L)))/(3.9*(v_ti**2)*tau_ii))
        #  ion-electron collision time (s):
        tau_ie = 1.5E-7*A*(Z**-2)*((Te)**(3/2))*(1/n)
        # tau_ie=1e-4;

        # Normalised ion-electron collision time:
        tau_ie_ = tau_ie/tau_n0
        # Normalised electron energy loss time:
        tau_Ee_ = ((3/5)*((1/tau_n + 1/tau_chi_e)**-1))/tau_n0
        # Normalised ion energy loss time:
        tau_Ei_ = ((3/5)*((1/tau_n + 1/tau_chi_i)**-1))/tau_n0
        yprime = np.array([-(n_*tau_n0)/tau_n,
                           -(Ee_/tau_Ee_)+((Ei_-Ee_)/tau_ie_),
                           -(Ei_/tau_Ei_)-((Ei_-Ee_)/tau_ie_)])
        return yprime

    def update_r(self, r0, tstep, v_elm, tau_n0):
        '''Function to update the position given the v_elm vector passed in
         by the user, the current time and the current position.
        '''
        # Get the current (initial) velocity and acceleration:
        if not np.shape(v_elm):
            v0 = v_elm
            acn = 0
        else:
            v0 = interp1d(v_elm[:, 0],
                          v_elm[:, 1])
            v0 = v0(r0)
            temp_1 = v_elm[:-1, 0]+0.5*np.diff(v_elm[:, 0]) # radial positions (halfway between specified points)
            temp_2 = np.diff(v_elm[:, 1])/np.diff(v_elm[:, 0]) # Delta v / Delta r , accelleration
            acn = interp1d(temp_1,
                           temp_2, kind='linear', fill_value='extrapolate') # accelleration function, requires a linear extrapolation because the first value of temp 1 is above zero, and then we try and intermp1 for the value at zero.
            acn = acn(r0) # accelleration found at position r0
        r_up = r0+v0*(tstep*tau_n0)+0.5*acn*(tstep*tau_n0)**2 # SUVAT equation for finding displacement in terms of initial velocity and accelleration
        return r_up

    def flx_lim(self, M, xi, v_ti, v_te):
        '''% Function to calculate the flux limiting factors alpha_e and
        alpha_i as a % function of M,xi,v_ti and v_te. Again, modifiable
        according to the regime. Set as in section 5.2.
        '''
        ae = 5*xi*(v_ti/v_te)
        ai = 3.5*xi
        return [ae, ai]

    def get_L(self, r):
        '''Function to calculate the connection length at a given time by
        interpolating the user-specified 2-column matrix, if given.
        '''
        if not np.shape(self.L):
            con = self.L
        else:
            # con = interp1d(self.L[:, 0], self.L[:, 1],
                               # fill_value='extrapolate')
            # con = con(r)
            con = self.interp1d_conLen(r)
        return con

    def ELM_loss_profile(self, v_elm, n0, Te0, Ti0, r, L):
        # Declare variables ###########################################################
        # Declare global variables so we don't have to keep passing them to f
        # global Z  # Proton number
        # global A  # Nucleon number
        # global tau_n0  # Initial parallel particle loss time
        # global Ee0  # Initial electron energy of the filament (MJ)
        # global Ei0  # Initial ion energy of the filament (MJ)
        # global n0  # Initial density of the filament (10^20 m^-3)
        # global MP  # Mass of a proton (kg)
        # global ME  # Mass of an electron (kg)
        # global QE  # Charge on an electron (Coulombs)
        v_elm = v_elm
        n0 = n0
        Te0 = Te0
        Ti0 = Ti0
        r_out = r
        n0 = n0
        Z = 1  # Deuterium plasma
        A = 2.5  # Deuterium - Tritium plasma
        ACC = 0.0005  # Desired accuracy
        MP = 1.673E-27  # Mass of a proton (kg)
        ME = 9.108E-31  # Mass of an electron (kg)
        QE = 1.602E-19  # Charge on an electron (Coulombs)

        out_len = len(r_out)
        n = UNN(out_len)
        Te = UNN(out_len)
        Ti = UNN(out_len)
        Ee = UNN(out_len)
        Ei = UNN(out_len)
        W = UNN(out_len)
        t_out = UNN(out_len)
        M = np.zeros(out_len)

        # Calculate the initial density loss time in SI units (Eq 5.3)
        xi0 = (1+Z/(Ti0/Te0))**0.5
        v_ti0 = ((Ti0*QE)/(MP*A))**0.5  # Initial ion thermal velocity
        tau_n0 = self.get_L(0)/(self.get_M(0)*xi0*v_ti0)
        Ee0 = (3/2)*n0*1E20*Te0*QE  # Initial average energy of the
                                    # filament in MJ
        Ei0 = (3/2)*n0*1E20*Ti0*QE  # Initial average energy of the
                                    # filament in MJ

        # Set the initial time step according to the maximum v_elm and the
        # minimum distance between successive values of r_out
        if len(r_out) > 1:
            if not np.shape(v_elm):
                tstep = (1/tau_n0)*np.min(np.diff(r_out))/np.max(v_elm)
            else:
                tstep = (1/tau_n0)*np.min(np.diff(r_out))/np.max(v_elm[:, 1])
        else:
            if not np.shape(v_elm):
                tstep = (1/tau_n0)*r_out/np.max(v_elm)
            else:
                tstep = (1/tau_n0)*r_out/np.max(v_elm[:, 1])

        # Dimensionless initial conditions in vector y (y=[n, Ee, Ei])
        y = np.array([1, 1, 1])

        # Start the time stepping. Time step until we have reached the last
        # position specified in r_out
        r = 0  # Starting position of the ELM filament
        t = 0  # Start at time zero

        # Prepare constants for the derivative function

        derivArgs = [Ee0, Ei0, QE, n0, A, MP, ME, Z, tau_n0]

        for r_out_i in range(len(r_out)):
            # Make a step (using RK4) that is either of acceptable length such
            # that the update is within specified accuracy or, if this length
            # oversteps a desired output position, step exactly to that output
            # position:
            rscout = r  # Scouts ahead until it goes past an output position
            tscout = t
            yscout = y
            while rscout < r_out[r_out_i]:
                # Update the current values of n', Ee' and Ei' using the RK4
                # method:
                r = rscout
                y = yscout
                t = tscout
                y1 = self.RK4(r, t, y, 2*tstep, self.f, derivArgs)  # One step of 2h
                # Two steps of h:
                y2 = self.RK4(r, t, y, tstep, self.f, derivArgs)
                r_temp = self.update_r(r, tstep, v_elm, tau_n0)
                y2 = self.RK4(r_temp, t, y2, tstep, self.f, derivArgs)
                delta = np.min(np.abs(y2-y1))  # Measure of accuracy
                # Update tstep according to eqn 16.2.7 of "Numerical Recipes in
                # C: the Art of Scientific Computing", chapter 16:
                if abs(delta) > ACC:
                    tstep = tstep*(ACC/delta)**0.2
                    continue
                else:
                    rscout = self.update_r(r_temp, tstep, v_elm, tau_n0)
                    tscout = tscout+2*tstep
                    yscout = y2
                    tstep = tstep*(ACC/delta)**0.2
            # If we're not exactly there now, update to the exact position at
            # which we want to calculate y:
            if rscout != r_out[r_out_i]:
                if not np.shape(v_elm):
                    tstep = (1/tau_n0)*(r_out[r_out_i]-r)/v_elm
                else:
                    tt = interp1d(v_elm[:, 0],
                                      v_elm[:, 1])
                    tt = tt(r)
                    tstep = (1/tau_n0)*(r_out[r_out_i]-r)/tt

                y = self.RK4(r, t, y, tstep, self.f, derivArgs)
                t += tstep
                r = self.update_r(r, tstep, v_elm, tau_n0)

            else:
                r = rscout
                y = yscout
                t = tscout

            n.UNNORM[r_out_i] = y[0]*n0
            n.NORM[r_out_i] = y[0]
            Ee.UNNORM[r_out_i] = y[1]*Ee0
            Ee.NORM[r_out_i] = y[1]
            Ei.UNNORM[r_out_i] = y[2]*Ei0
            Ei.NORM[r_out_i] = y[2]
            Te.UNNORM[r_out_i] = (2/3)*(Ee.UNNORM[r_out_i] /
                                        (n.UNNORM[r_out_i]*1E20*QE))
            Te.NORM[r_out_i] = Te.UNNORM[r_out_i]/Te0
            Ti.UNNORM[r_out_i] = (2/3)*(Ei.UNNORM[r_out_i] /
                                        (n.UNNORM[r_out_i]*1E20*QE))
            Ti.NORM[r_out_i] = Ti.UNNORM[r_out_i]/Ti0
            W.UNNORM[r_out_i] = (Ee.UNNORM[r_out_i]+Ei.UNNORM[r_out_i])
            W.NORM[r_out_i] = (y[1]+y[2])/2
            t_out.UNNORM[r_out_i] = t*tau_n0
            t_out.NORM[r_out_i] = t
            M[r_out_i] = self.get_M(t)

        return [n, Te, Ti, Ee, Ei, W, t_out, M]

    def calculate_ELM_profile(self, W_harm_med, r, WELM0, fELM,
                              nfil, poloWidth, aveFac, t_out, v_elm,
                              Etot, sigmar, sigmaz, lambdaW,
                              dtout2=1e-6, drsep=[0, 0.02, 0.04, 0.06,
                                                  0.08, 0.1],
                              rg2f=2, npts=500, rlim=0.3):

        WELM = W_harm_med.NORM  # fraction of W0
        t_out = t_out.UNNORM # time
        drsep1 = t_out*v_elm # time converted to sep1 distance

        # (1) Wfil as a function of time
        t_out2 = np.arange(np.min(t_out), np.max(t_out), dtout2)
        f_WELM2_t = interp1d(t_out, WELM, 'cubic') #making the time
                                                   #trace of Wfil
                                                   #smoother
        WELM2 = f_WELM2_t(t_out2)
        ntout2 = len(t_out2)

        # (2) q|| instant. as a function of time
        dWdt=[]
        tmid=[]

        for ii in range(ntout2-1):
            tmid.append((t_out2[ii]+t_out2[ii+1])/2)
            dWdt.append( -(WELM2[ii+1]-WELM2[ii])/(t_out2[ii+1]-t_out2[ii]))

        qcenter = []
        for ii in range(ntout2-1):
            qcenter.append(0.5*(Etot/nfil)*dWdt[ii]/(2*np.pi*sigmar*sigmaz))

        rr = np.linspace(0, rlim, npts)

        Epar = np.zeros(len(rr))   # numerical integration with
                                   # spatially varying lambda
        EparAn = np.zeros(len(rr)) # analytical approximation for
                                   # constant lambda
        for jj in range(len(rr)):
            qpar=[];
            for ii in range(ntout2-2):
                qpar.append(0.5*(Etot/nfil)*dWdt[ii] \
                            *np.exp(-(rr[jj] -
                                      tmid[ii]*v_elm)**2 /
                                    (2*sigmar**2))/(2*np.pi*sigmar*sigmaz))

            for ii in range(ntout2-3):
                Epar[jj] = Epar[jj] + (qpar[ii]+qpar[ii++1])/2 * \
                    (t_out2[ii+1]-t_out2[ii])

            EparAn[jj]=0.25*Etot/(nfil*lambdaW*np.sqrt(2*np.pi)*sigmaz) \
            * np.exp(sigmar**2/(2*lambdaW**2)) * \
            np.exp(-rr[jj]/lambdaW) * \
            (1 + math.erf((rr[jj]/sigmar -
                           sigmar/lambdaW)/np.sqrt(2)))

        qpar = Epar*fELM/rg2f/1e6

        # (7) characteristic e-folding length of WELM
        LL=[]
        rr2=[]
        for ii in range(len(drsep1)-1):
            rr2.append((drsep1[ii+1]+drsep1[ii])/2)
            LL.append(-(((WELM[ii+1]-WELM[ii])/(drsep1[ii+1]  \
                        - drsep1[ii]))/ \
                        ((WELM[ii+1]+WELM[ii])/2))**-1)

        return [rr, qpar, LL, WELM2, qcenter, Epar, EparAn, t_out2]


def get_elm_data(conlen_graph_data, generate_graphics=False, output_name="",
    Rb=1, Btot=1, Bpm=1, r_break=0.025, drsep=1.0):
    """Performs the ELM PLM calculation.
    """
    r_out = conlen_graph_data[:, 0]
    r_out[0] = 0 # Set the first drsep to 0 if it isn't...
    Lcon_data = conlen_graph_data[:, 3]

    fELM = 33 # ELM frequency [Hz]
    WELM0 = 0.6  # Energy lost per ELM [MJ]
    nfil = 10  # number of filaments
    poloWidth = 0.6  # poloidal width of filament [m]
    aveFac = 1.7  # averaging factor
    Aion = 2.5  # 2.5 for DT plasma
    n0 = 0.375  # filament peak density at R=Rsep [10^20 m^-3]
    Te0 = 2500  # filament peak electron temperature at R=Rsep [eV]
    Ti0 = Te0  # filament peak ion temperature at R=Rsep [eV]
    dr = 0  # difference of r from the LCFS
    v_elm = 500  # radial propagation speed [m/s]
    Etot = 0.6e6  # prescribed total energy released per ELM, dWELM [W]
    sigmar = 3e-2  # filament radial width [m]
    sigmaz = 0.45  # filament vertical height [m]
    rg2f = 2  # gap-to-filament ratio
    nfil = 10  # number of ELM filaments
    lambdaW = 0.0213  # approximate radial e-folding length of the ELM
                      # energy (fit for drsep1 = 0-10 cm for the QDN
                      # equilibrium)
    log.info("Calculating ELM profile")
    elm_loss = ELM_loss(v_elm, n0, Te0, Ti0, Aion,
                        fELM, nfil, poloWidth, aveFac, Lcon_data, r_out)
    [n, Te, Ti, Ee, Ei,
     W, t_out, M] = elm_loss.ELM_loss_profile(v_elm, n0, Te0, Ti0, r_out,
                                              Lcon_data)
    W_harm_med = W

    data_jnme = elm_loss.calculate_ELM_profile(W_harm_med,
                                               r_out,
                                               WELM0, fELM,
                                               nfil,
                                               poloWidth,
                                               aveFac, t_out,
                                               v_elm, Etot,
                                               sigmar,
                                               sigmaz,
                                               lambdaW,
                                               npts=len(r_out),
                                               rlim=max(r_out))


    x = np.empty((len(data_jnme[0]), 2), np.float64)
    x[:, 0] = data_jnme[0]
    x[:, 1] = data_jnme[1]
    log.info("Finished calculating ELM profile")

    if not generate_graphics:
        # Just return the value
        return x

    import matplotlib.pyplot as plt

    log.info("Generating ELM graphics")

    fig, axs = plt.subplots(1, 4, figsize=(12, 6), tight_layout=True, dpi=100)

    axs[0].set_title('a) $L_{conn}[m]$')
    axs[1].set_title('b) $W_{fil}/W_{0fil}$')
    axs[2].set_title('c) $\lambda [mm]$')  # $[(dW/dr)/W]^{-1}[mm]$
    axs[3].set_title(r'a) $q_{\parallel, ELM} + q_{\parallel, inter-ELM}[\frac{MW}{m^2}]$')

    XLABEL = '$r-r_{sep}[mm]$'

    axs[0].plot(r_out * 1e3, Lcon_data)
    axs[0].grid(True)
    axs[0].set_xlabel(XLABEL)
    axs[0].set_xlim((-2.5, 200.0))

    axs[1].grid(True)
    axs[1].set_xlim((0, 200.0))
    axs[1].set_xlabel(XLABEL)
    axs[1].semilogy(r_out[:-1]*1000, W_harm_med.NORM[:-1], 'b',
                  label='Steady state')

    axs[2].semilogy(np.array(data_jnme[0][:-1]) * 1000,
                    np.array(data_jnme[2]) * 1000, 'b',
                    label='Steady state')
    axs[2].grid(True)
    axs[2].set_xlabel(XLABEL)
    axs[2].set_xlim((-2.5, 200.0))

    # Get interELM in MW/m^2 !
    q_interELM_break = inter_ELM(data_jnme[0], Rb, Btot, Bpm, Rb=r_break) * 1e-6

    axs[3].semilogy(data_jnme[0] * 1000, data_jnme[1] + q_interELM_break, 'b',
                    label='ELM PLM + inter-ELM')
    axs[3].semilogy(data_jnme[0] * 1000, data_jnme[1], 'y--',
                    label='ELM PLM')
    axs[3].semilogy(data_jnme[0] * 1000, q_interELM_break, 'g--',
                    label='inter-ELM')

    axs[3].set_xlim((0, 200.0))
    axs[3].set_ylim((0.001, 1000.0))

    # Steady-State profile
    DRSEP = drsep * 1e-3 # We convert it back to meters
    r = np.linspace(DRSEP, 0.3, 500)
    # 8.29 - 8.20

    q_elm = 5 * np.exp(-(r - DRSEP) / 0.09)
    q_inter = 3 * np.exp(-(r - DRSEP) / 0.17)


    axs[3].semilogy(r*1000, q_elm + q_inter, 'k', label='Conservative model')
    # axs[3].semilogy(r*1000 + 90, q_elm, 'k', label='Conservative model')
    # plt.tick_params(axis='y', which='minor')
    # axs[3].yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    # axs[3].set_ylabel(r'$q_{\parallel} \frac{MW}{m^2}$'))
    axs[3].set_xlabel(XLABEL)
    axs[3].grid(True, which='both')
    axs[3].vlines(1e3 * DRSEP, ymin=0, ymax=200, colors='r')
    axs[3].text(1e3 * DRSEP, 150, r"$2^{nd} sep$", rotation=90)
    axs[3].legend()

    fig.savefig(f"{output_name}.pdf")

    # Now just a figure with the heat loads profiles.

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.semilogy(data_jnme[0] * 1000, data_jnme[1] + q_interELM_break, 'b',
                label='ELM PLM + inter-ELM')
    ax.semilogy(data_jnme[0] * 1000, data_jnme[1], 'y--',
                label='ELM PLM')
    ax.semilogy(data_jnme[0] * 1000, q_interELM_break, 'g--',
                label='inter-ELM')
    ax.semilogy(r*1000, q_elm + q_inter, 'k', label='Conservative model')
    ax.set_xlabel(XLABEL)
    ax.grid(True, which='both')
    ax.vlines(1e3 * DRSEP, ymin=0, ymax=200, colors='r')
    ax.text(1e3 * DRSEP, 150, r"$2^{nd} sep$", rotation=90)
    ax.set_title(r'a) $q_{\parallel, ELM} + q_{\parallel, inter-ELM}[\frac{MW}{m^2}]$')
    ax.legend()
    ax.set_xlim((0, 200.0))
    ax.set_ylim((0.001, 1000.0))

    f.savefig(f"{output_name}_qpar.pdf")

    return x