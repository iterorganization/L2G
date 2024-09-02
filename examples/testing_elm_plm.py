import l2g.hlm.elm_plm
import l2g.hlm.steady_state
import l2g.equil
import l2g.external.equilibrium_analysis

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


    eq_file = "shot_135011_run_7_t_399.927598s.eqdsk"

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


    obj = l2g.hlm.elm_plm.ELM_PLM(drsep, conlen_down)
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