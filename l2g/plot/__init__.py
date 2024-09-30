#!/usr/bin/env python3

from l2g.plot._plot_psi import plot_psi

def plotPolyLine(ax, contour, label='', color='k'):
    x,y = contour

    labelFlag = 1
    N = len(x)

    lines = []
    for path in range(N):
        pX = x[path]
        pY = y[path]
        if labelFlag:
            lines += ax.plot(pX, pY, color=color, label=label)
            labelFlag = 0
        else:
            lines += ax.plot(pX, pY, color=color)
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
    from matplotlib.animation import FuncAnimation, PillowWriter
    from functools import partial
    import sys # For printing to stdout

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
