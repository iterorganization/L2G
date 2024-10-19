#!/usr/bin/env python3

from l2g.plot._marching_squares import Marching
from l2g.plot._polygon import Polygon

def plotPolyLine(ax, contour, label='', color='k', linewidth=1.0):
    x,y = contour

    labelFlag = 1
    N = len(x)

    lines = []
    for path in range(N):
        pX = x[path]
        pY = y[path]
        if labelFlag:
            lines += ax.plot(pX, pY, color=color, label=label,
                             linewidth=linewidth)
            labelFlag = 0
        else:
            lines += ax.plot(pX, pY, color=color, linewidth=linewidth)
    return lines

def plotContour(equil: 'Equilibrium', ax, psi_boundary=None, levels_core=None, levels_vacc=None):
    from scipy.interpolate import interp2d

    I_Psi = interp2d(equil.grid_r, equil.grid_z, equil.psi, kind='cubic')

    ax.contour(equil.grid_r, equil.grid_z,
               I_Psi(equil.grid_r, equil.grid_z, 0, 0), 35, cmap="jet")

def plotLimiter(equil: 'Equilibrium', ax):
    ax.plot(equil.wall_contour_r, equil.wall_contour_z, 'r-')

def plotQuiever(equil: 'Equilibrium', ax, res_width=20, res_height=40):
    """Quiever plot.

    res_width is the number of arrows by width
    res_height is the number of arrows by height
    scale is some matplotlib non-linear non-deterministic scale for arrows...
    """
    # Quiever plot of the poloidal components

    if equil.grid_dim_r > res_width:
        nwStep = equil.grid_dim_r // res_width
    else:
        nwStep = 1
    if equil.grid_dim_z > res_height:
        nhStep = equil.grid_dim_z // res_height
    else:
        nhStep = 1

    R = equil.grid_r[::nwStep]
    Z = equil.grid_z[::nhStep]

    from scipy.interpolate import interp2d
    I_Psi = interp2d(equil.grid_r, equil.grid_z, equil.psi, kind='cubic')
    BR = -I_Psi(R, Z, 0, 1) / R
    BZ = I_Psi(R, Z, 1, 0) / R

    ax.quiver(R, Z, BR, BZ, pivot='tip')


def plotEqdskLcfs(eq: 'EQ', ax, color1='g', color2='b', label=''):
    # eq.evaluate()

    lcfs1 = eq.getContourPaths(flux=eq.psiLCFS)
    plotPolyLine(ax, lcfs1, label, color=color1)
    if eq.type_ == 'div':
        lcfs2 = eq.getContourPaths(flux=eq.psiLCFS2)
        plotPolyLine(ax, lcfs2, '', color=color2)

def plotOPoint(equil: 'Equilibrium', ax, pattern='ro'):
    ax.plot(equil.mag_axis_r, equil.mag_axis_z, pattern)

def plotMidplane(equil: 'Equilibrium', ax, pattern='r-'):
    zmaxis = equil.mag_axis_z
    rleft = equil.grid_r[0]
    rend = equil.grid_r[-1]
    ax.plot([rleft - 1.5, rend + 1.5], [zmaxis, zmaxis], pattern)

def createLcfsMovie(equilibriums, labels, legend_name='Legend',
    legend_pos='upper right', movieName='movie.gif', dpi=90, bitrate=1000,
    type='gif'):
    """Creates a movie out of EQDSKs. It takes a list of equilibriums, takes
    the LCFS or boundary points and plot it.
    """
    if not isinstance(equilibriums, list):
        equilibriums = [equilibriums]

    # if not labels:
    #     labels = [_.getName() for _ in eqdsk]
    from l2g.equil import EQ
    import matplotlib.pyplot as plt
    import matplotlib.figure
    import matplotlib.axes
    from matplotlib.animation import FuncAnimation, PillowWriter
    from functools import partial
    import sys # For printing to stdout

    fig: matplotlib.figure.Figure
    ax: matplotlib.axes.Axes

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    ax.grid(True)
    ax.axis('equal')

    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')

    plotLimiter(equilibriums[0], ax)

    # Plot first LCFS
    eq = EQ()
    eq.setEquilibrium(equilibriums[0])
    eq.evaluate()
    path = eq.getContourPaths(flux=eq.psiLCFS)
    lines = plotPolyLine(ax, path, label='test')
    while lines:
        # The pop(0) element gets us the ax.plot object. And to completely
        # remove it you also have to call remove()
        lines.pop(0).remove()

    # xlimits = ax.get_xlim()
    # ylimits = ax.get_ylim()
    fig.canvas.draw()
    N = len(equilibriums)

    sys.stdout.write(f"\rProcessing movie: {int(100 * 0.0/ N)}%")
    sys.stdout.flush()
    if type == 'gif':
        def update(i, lines):
            sys.stdout.write(f"\rProcessing movie: {100 * i/ N:.2f}%")
            sys.stdout.flush()
            while lines:
                lines.pop(0).remove()

            eq.setEquilibrium(equilibriums[i])
            eq.evaluate()
            path = eq.getContourPaths(flux=eq.psiLCFS)
            lines += plotPolyLine(ax, path, label=labels[i], color='b')

            if eq.psiLCFS2:
                path = eq.getContourPaths(flux=eq.psiLCFS2)
                lines += plotPolyLine(ax, path, color='g')
            ax.legend(title=legend_name, loc=legend_pos)

        if not movieName.lower().endswith('.gif'):
            movieName += '.gif'
        anim = FuncAnimation(fig, partial(update, lines=lines), frames=N, interval=250)
        anim.save(movieName, writer=PillowWriter(bitrate=bitrate))
    else:
        # Plot the LCFSs then save them in an output
        for i in range(N):
            sys.stdout.write(f"\rProcessing movie: {100 * i/ N:.2f}%")
            sys.stdout.flush()
            while lines:
                lines.pop(0).remove()
            eq.setEquilibrium(equilibriums[i])
            eq.evaluate()
            path = eq.getContourPaths(flux=eq.psiLCFS)
            lines += plotPolyLine(ax, path, label=labels[i], color='b')

            if eq.psiLCFS2:
                path = eq.getContourPaths(flux=eq.psiLCFS2)
                lines += plotPolyLine(ax, path, color='g')
            ax.legend(title=legend_name)

                # Save the image as png.
            if not movieName.lower().endswith('.png'):
                out = movieName + f'_{i}.png'
            else:
                out = movieName[:-4] + f'_{i}.png'
            fig.savefig(f'{out}')

def plot_psi_to_mpl_ax(ax, eq: 'l2g.equil.EQ'):
    """Plot the poloidal magnetic flux of an equilibrium to a matplotlib axes.

    Arguments:
        ax: matplotlib axis or subplot
        eq: Equilibrium with diagnostic
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np

    eq.evaluate()
    psi = eq._eq.psi
    grid_r = eq._eq.grid_r
    grid_z = eq._eq.grid_z


    psi_axis = eq._eq.psi_axis
    psi_boundary = eq.psiLCFS
    psi_2nd_boundary = eq.psiLCFS2 # In case of limiter this is None
    equ_type = eq.type_
    if equ_type == "div" and psi_2nd_boundary is None:
        # Diverted with single X point
        equ_type = "lim"

    psi_outside = eq._psi_spline.ev(grid_r[-1], eq._eq.mag_axis_z)

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    base_cmap_core = plt.get_cmap('jet_r')
    base_cmap_drsep = plt.get_cmap('jet_r')
    base_cmap_vacc = plt.get_cmap('jet_r')

    # Create intervals
    cmap_vacc = truncate_colormap(base_cmap_vacc, minval=0.8, maxval=1.0, n=10)
    if equ_type == "div":
        # VACUUM

        if psi_2nd_boundary > psi_outside:
            band = np.linspace(psi_outside, psi_2nd_boundary, 10)
        else:
            band = np.linspace(psi_2nd_boundary, psi_outside, 10)

        ax.contour(grid_r, grid_z, psi, levels=band, cmap=cmap_vacc, alpha=0.3)
        ax.contourf(grid_r, grid_z, psi, levels=band, cmap=cmap_vacc,
                    alpha=0.2)
        # DRSEP
        if psi_boundary > psi_2nd_boundary:
            band = np.linspace(psi_2nd_boundary, psi_boundary, 10)
        else:
            band = np.linspace(psi_boundary, psi_2nd_boundary, 10)

        cmap_drsep = truncate_colormap(base_cmap_drsep, minval=0.5, maxval=0.8, n=10)
        ax.contour(grid_r, grid_z, psi, levels=band, cmap=cmap_drsep,
                   alpha=0.5)
        ax.contourf(grid_r, grid_z, psi, levels=band, cmap=cmap_drsep,
                    alpha=0.2)
    else:
        # VACUUM
        if psi_boundary > psi_outside:
            band = np.linspace(psi_outside, psi_boundary, 10)
        else:
            band = np.linspace(psi_boundary, psi_outside)
        ax.contour(grid_r, grid_z, psi, levels=band, cmap=cmap_vacc, alpha=0.3)
        ax.contourf(grid_r, grid_z, psi, levels=band, cmap=cmap_vacc,
                    alpha=0.2)


    # CORE
    cmap_core = truncate_colormap(base_cmap_core, minval=0, maxval=0.3, n=25)
    if psi_axis > psi_boundary:
        band = np.linspace(psi_boundary, psi_axis, 25)
    else:
        band = np.linspace(psi_axis, psi_boundary, 25)

    ax.contour(grid_r, grid_z, psi, levels=band, cmap=cmap_core)
    ax.contourf(grid_r, grid_z, psi, levels=band, cmap=cmap_core, alpha=0.2)

    lcfs1 = eq.getContourPaths(flux=eq.psiLCFS)
    if equ_type == 'div':
        plotPolyLine(ax, lcfs1, label=r"$1^{st}$", color='g', linewidth=1.5)
        lcfs2 = eq.getContourPaths(flux=eq.psiLCFS2)
        plotPolyLine(ax, lcfs2, label=r"$2^{nd}$", color='b', linewidth=1.5)
        # 4cm contour.
        Rbdry, Zcen, _, _ = eq.getOWL_midplane()

        psi_dist4cm = eq._psi_spline.ev(Rbdry + 0.04, Zcen)
        con_dist4cm = eq.getContourPaths(flux=psi_dist4cm)
        plotPolyLine(ax, con_dist4cm, label=r"40 mm flux surface", color='tab:orange', linewidth=1.5)
    else:
        plotPolyLine(ax, lcfs1, label=r"LCFS", color='g', linewidth=1.5)

    # Plot limiter silhouette
    ax.plot(eq._eq.wall_contour_r, eq._eq.wall_contour_z, 'r-', linewidth=2.5)
    ax.legend()