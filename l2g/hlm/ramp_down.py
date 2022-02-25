
def decay_length_L_mode_diverted(a: float, R: float, Ip: float, Area: float):
    r"""Calculates the decay length of a L-mode diverted equilibrium, by the
    use of the most suitable scaling law from Horacek, et.al, 2020.

    For convenience the data stored in a IMAS IDS is in the following:

     - a = equilibirum.time_slice[:].boundary.minor_radius
     - R = equilibrium.time_slice[:].global_quantities.magnetic_axis.r
     - Ip = equilibirum.time_slice[:].global_quantities.ip
     - Area = equilibrium.time_slice[:].global_quantities.area

    .. math::

       \lambda_q = 4350 (\frac{a}{R})^{1.09} (\frac{Ip}{Area})^{-0.43}


    Returns:
        lambda_q (float): Decay length. In millimeters.

    Arguments:
        a (float): minor radius or plasma width at the midplane. In meters.
        R (float): Magnetic axis major radius. In meters.
        Ip (float): Plasma current. In MA.
        Area (float): Plasma cross section. In m^2.

    """

    lambda_q = 4350 * (a / R)**(1.09) * (Ip / Area)**(-0.43)

    return lambda_q