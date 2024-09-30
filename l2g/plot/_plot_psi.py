"""This file contains the function to plot the magnetic flux map of an
equilibrium.
"""
import os

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.figure
import numpy as np

def plot_psi_to_mpl_ax(ax, eq: 'l2g.equil.EQ'):
    """Plot the poloidal magnetic flux of an equilibrium to a matplotlib axes.

    Arguments:
        ax: matplotlib axis or subplot
        eq: Equilibrium with diagnostic
    """
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

    def plotPolyLine(ax, contour, label='', color='k', linewidth=1):
        x,y = contour

        labelFlag = 1
        N = len(x)

        lines = []
        for path in range(N):
            pX = x[path]
            pY = y[path]
            if labelFlag:
                lines += ax.plot(pX, pY, color=color, label=label, linewidth=linewidth)
                labelFlag = 0
            else:
                lines += ax.plot(pX, pY, color=color, linewidth=linewidth)
        return lines

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


def plot_psi(eq: 'l2g.equil.EQ', output_path: str = "",
             plot_midplane: bool = False,
             save_files: bool = True) -> matplotlib.figure.Figure:
    """Plot the poloidal magnetic flux of an equilibrium.

    If for output only the directory is specified, it generates two graphics in
    the directory, named psi.png and psi.pdf, otherwise it will append the
    _psi and generate a pdf and png file.

    Arguments:
        eq (l2g.equil.EQ'): Equilibrium with diagnostic
        output_path (str): Output path.
    """

    figure = plt.figure(figsize=plt.figaspect(3/2))
    ax = figure.add_subplot(111)
    ax.axis("equal")
    ax.set_xlabel("R[m]")
    ax.set_ylabel("Z[m]")
    plot_psi_to_mpl_ax(ax, eq)

    figure.tight_layout()

    if save_files:
        if os.path.isdir(output_path):
            pdf_name = os.path.join(output_path, "psi.pdf")
            png_name = os.path.join(output_path, "psi.png")
        else:
            pdf_name = f"{output_path}_psi.pdf"
            png_name = f"{output_path}_psi.png"

        figure.savefig(pdf_name)
        figure.savefig(png_name)
    return figure