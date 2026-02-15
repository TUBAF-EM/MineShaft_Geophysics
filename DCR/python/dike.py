import numpy as np


def dike2(RHO1, RHO2, THICK, YC, YP, XC=0, XP=0):
    """
    Pole-pole model response for vertical dike.

    Parameters
    ----------
    RHO1 : float
        Resistivity of left vertical layer
    RHO2 : float
        Resistivity of the middle vertical layer (dike)
    THICK : float
        Thickness of the middle layer (dike)
    XC, YC : float
        Coordinates of current electrode
    XP, YP : float
        Coordinates of potential electrode

    Returns
    -------
    iterable
        Voltage measured at the potential electrode
    """

    MAX = 100
    PI2 = 2 * np.pi

    X = XP - XC
    X2 = X * X
    THICK2 = THICK + THICK
    AK = (RHO2 - RHO1) / (RHO2 + RHO1)
    AK2 = AK * AK
    CC = 1.0 + AK
    DD = 1.0 - AK2
    AA = -AK * DD
    BB = -AK * CC

    # Same position → zero potential
    if np.all(np.isclose(XC, XP) & np.isclose(YC, YP)):
        return 0.0

    # ------------------------------------------------------------------
    # CASE 1: YC <= 0
    # ------------------------------------------------------------------
    elif YC <= 0.0:
        S = abs(YC)
        S2 = S + S

        if YP <= 0.0:
            A = YP - YC
            Y = (MAX + 1) * THICK2 + S2 - A
            SUM = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = M * THICK2 + S2 - A
                SUM = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM

            V = (
                1.0 / np.sqrt(X2 + A * A)
                + AK / np.sqrt(X2 + (S2 - A) ** 2)
                + AA * SUM
            )
            return RHO1 / PI2 * V

        elif YP <= THICK:
            A = abs(YP) + abs(YC)

            Y = (MAX + 1) * THICK2 + S2 - A
            SUM1 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = MAX * THICK2 + A
            SUM2 = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = M * THICK2 + S2 - A
                SUM1 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM1

                Y = (M - 1) * THICK2 + A
                SUM2 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM2

            V = BB * SUM1 + CC * SUM2
            return RHO1 / PI2 * V

        else:
            A = abs(YP) + abs(YC)
            Y = MAX * THICK2 + A
            SUM = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = (M - 1) * THICK2 + A
                SUM = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM

            V = DD * SUM
            return RHO1 / PI2 * V

    # ------------------------------------------------------------------
    # CASE 2: 0 < YC <= THICK
    # ------------------------------------------------------------------
    elif (YC > 0.0) and (YC <= THICK):
        S = YC
        S2 = S + S

        if YP <= 0.0:
            A = abs(YP) + YC

            Y = MAX * THICK2 + A
            SUM1 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = (MAX + 1) * THICK2 - S2 + A
            SUM2 = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = (M - 1) * THICK2 + A
                SUM1 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM1

                Y = M * THICK2 - S2 + A
                SUM2 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM2

            V = (1.0 + AK) * (SUM1 - AK * SUM2)
            return RHO1 / PI2 * V

        elif YP <= THICK:
            A = YP - YC

            Y = (MAX + 1) * THICK2 - A
            SUM1 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = (MAX + 1) * THICK2 + A
            SUM2 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = MAX * THICK2 + S2 + A
            SUM3 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = (MAX + 1) * THICK2 - S2 - A
            SUM4 = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = M * THICK2 - A
                SUM1 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM1

                Y = M * THICK2 + A
                SUM2 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM2

                Y = (M - 1) * THICK2 + S2 + A
                SUM3 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM3

                Y = M * THICK2 - S2 - A
                SUM4 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM4

            V = (
                1.0 / np.sqrt(X2 + A * A)
                + AK2 * (SUM1 + SUM2)
                - AK * (SUM3 + SUM4)
            )
            return RHO2 / PI2 * V

        else:
            A = YP - YC

            Y = MAX * THICK2 + A
            SUM1 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = MAX * THICK2 + S2 + A
            SUM2 = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = (M - 1) * THICK2 + A
                SUM1 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM1

                Y = (M - 1) * THICK2 + S2 + A
                SUM2 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM2

            V = (1.0 + AK) * (SUM1 - AK * SUM2)
            return RHO1 / PI2 * V

    # ------------------------------------------------------------------
    # CASE 3: YC > THICK
    # ------------------------------------------------------------------
    else:
        S = YC - THICK
        S2 = S + S

        if YP <= 0.0:
            A = abs(YP) + YC
            Y = MAX * THICK2 + A
            SUM = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = (M - 1) * THICK2 + A
                SUM = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM

            V = DD * SUM
            return RHO1 / PI2 * V

        elif YP <= THICK:

            A = abs(YP - YC)

            Y = MAX * THICK2 + A
            SUM1 = 1.0 / np.sqrt(X2 + Y * Y)

            Y = (MAX + 1) * THICK2 + S2 - A
            SUM2 = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = (M - 1) * THICK2 + A
                SUM1 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM1

                Y = M * THICK2 + S2 - A
                SUM2 = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM2

            V = CC * (SUM1 - AK * SUM2)
            return RHO1 / PI2 * V

        else:
            A = YP - YC

            Y = (MAX + 1) * THICK2 + S2 + A
            SUM = 1.0 / np.sqrt(X2 + Y * Y)

            for M in range(MAX, 0, -1):
                Y = M * THICK2 + S2 + A
                SUM = 1.0 / np.sqrt(X2 + Y * Y) + AK2 * SUM

            V = (1.0 / np.sqrt(X2 + A * A) +
                 AK / np.sqrt(X2 + (S2 + A) ** 2) +
                 AA * SUM)

            return RHO1 / PI2 * V

def getDike(rho1, rho2, x, thick, pos=0):
    """
    Model response for vertical dike for different electrode arrays.

    Parameters
    ----------
    rho1 : float
    rho2 : float
    x : array_like
        Measurement positions (1D array)
    thick : float
    pos: float
        Position of the dike (start)

    Returns
    -------
    rhosKN, rhosSB, rhosOK, rhosUK : numpy.ndarray
        Apparent resistivities for different arrays (KN, SL, OK, UK)
    """

    y = np.asarray(x, dtype=float) - pos
    na = len(y)

    rhosKN = np.zeros_like(y)
    rhosSB = np.zeros_like(y)
    rhosOK = np.zeros_like(y)
    rhosUK = np.zeros_like(y)

    anordnung = ['KN', 'SL', 'OK', 'UK']

    for cur_anordnung in anordnung:
        # electrode distances
        if cur_anordnung == 'KN':
            am = 0.1

        elif cur_anordnung in ['SL', 'OK', 'UK']:
            am = 0.4
            mn = 0.1

        # Iteration over measurement positions
        for ii in range(na):
            if cur_anordnung == 'KN':
                um = dike2(rho1, rho2, thick,
                           y[ii] + am/2,
                           y[ii] - am/2)

                rhosKN[ii] = um * 2*np.pi / (1.0/am)

            elif cur_anordnung == 'SL':
                uam = dike2(rho1, rho2, thick,
                            y[ii] + (am + mn/2),
                            y[ii] + mn/2)

                uan = dike2(rho1, rho2, thick,
                            y[ii] + (am + mn/2),
                            y[ii] - mn/2)

                ubm = dike2(rho1, rho2, thick,
                            y[ii] - (am + mn/2),
                            y[ii] - mn/2)

                ubn = dike2(rho1, rho2, thick,
                            y[ii] - (am + mn/2),
                            y[ii] + mn/2)

                factor = 2*np.pi / (1.0/am - 1.0/(am + mn))

                rhosSB[ii] = 0.5 * (
                    (uam - uan) * factor +
                    (ubm - ubn) * factor
                )

            elif cur_anordnung == 'OK':

                um = dike2(rho1, rho2, thick,
                           y[ii] + (am + mn/2),
                           y[ii] + mn/2)

                un = dike2(rho1, rho2, thick,
                           y[ii] + (am + mn/2),
                           y[ii] - mn/2)

                rhosOK[ii] = (um - un) * 2*np.pi / (1.0/am - 1.0/(am + mn))

            elif cur_anordnung == 'UK':

                um = dike2(rho1, rho2, thick,
                           y[ii] - (am + mn/2),
                           y[ii] - mn/2)

                un = dike2(rho1, rho2, thick,
                           y[ii] - (am + mn/2),
                           y[ii] + mn/2)

                rhosUK[ii] = (um - un) * 2*np.pi / (1.0/am - 1.0/(am + mn))

    return rhosKN, rhosSB, rhosOK, rhosUK
