def T_eval(Q_edgeBe, d):
    """Temperature evaluation by https://user.iter.org/default.aspx?uid=2ZPSMQ

    T = T_ref + ((Q_edgeBe - 0.1762) * d + 0.419)/(0.0265 + .0113d)
    Where:
        T_ref = 615 degrees
        Q_edge is heat load on Beryllium tile in MW/m^2
        d is depth in mm
    """

    return 615 + ((Q_edgeBe - 0.1762)*d + 0.419)/(0.0265 + 0.0113 * d)
"""
16x16 tile

P = 120W
T ~= 800C
"""

def P_eval(Q_edgeBe, d, tile=16):
    """
    Power evaluation by https://user.iter.org/default.aspx?uid=TFK7TJ for 12x12 and 16x16 mm tiles.

    """

    if tile == 16:
        power = d * 16 * Q_edgeBe
        maxP = 120
    elif tile == 12:
        power = d * 12 * Q_edgeBe
        maxP = 68
    else:
        print(f"Wrong tile: {tile}")
        return

    percentage = power / maxP

    print(f"Power: {power:.2f} W; {percentage*100:.1f}% of max power: {maxP}W")
    return power


def max_Q(d):
    """ Based on max proposed edge heat load value from https://user.iter.org/default.aspx?uid=TFK7TJ
    calculate what would be the maximum allowed heat load on certain depth.

    Q_edgeBe = P/d**2

    Where:
        P=68W - power in wats
        d - depth in mm
        Q_edgeBe - heat load in MW/m^2

    The maximum allowed power was evaluated using a 12x12 mm2 tile on which the
    top heat load is 680W. Then the top heat load shall be 10% of that. The rationale
    as quoted from the IDM document is:
    ```
    The rationale for this 10% criterion is that the added stress at the
    bonding interface remains small, compared to the existing interfacial
    stress that is already present because of the heat load passing through
    the structure
    ```

    Similarly the proposed criterion can be defined as line power in W/m, which nets
    to:
        k = 5.7 kW/m

    Resulting in the following equation:

    Q_edgeBe = k / d


    """

    # Check if P/d**2 is almost the same as 5.7 * d, where d is in mm.

    return 5.7 / d,  68 / (d*12), 7.5 / d, 120 / (d * 16)

def plot_MaxQ_ofD():
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.linspace(0.1, 50, 100)
    y = max_Q(x)
    plt.plot(x, y)

    plt.xlabel('Depth [mm]')
    plt.ylabel(r'$Q_{edge-Be} [\frac{MW}{m^2}]$')
    plt.yscale('log')
    plt.grid(True)
    plt.show()
