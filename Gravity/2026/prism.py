import numpy as np

def gz_prism(obs=(0.0, 0.0, 0.0), prism=(1.0, 2.0, 1.0, 2.0, 1.0, 2.0), rho=2700.0):
    """
    Calculate the vertical gravity component gz (after Nagy 1966)
    at point (x, y, z) for a homogeneous rectangular prism.

    Parameters
    ----------
    obs : tuple
        Observation point x, y, z (in m)
    prism : tuple
        Prism boundaries (x1, x2, y1, y2, z1, z2) in m
        (z1 and z2 relative to the same reference level as z)
    rho : float
        Density contrast (kg/m³)

    Returns
    -------
    gz : float
        Vertical gravity component (in mGal)
    """
    G = 6.67430e-11
    x, y, z = obs
    x1, x2, y1, y2, z1, z2 = prism
    gz = 0.0

    for i, xi in enumerate([x1 - x, x2 - x]):
        for j, yj in enumerate([y1 - y, y2 - y]):
            for k, zk in enumerate([z1 - z, z2 - z]):
                sgn = (-1)**(i + j + k)
                r = np.sqrt(xi**2 + yj**2 + zk**2)
                # if r == 0:
                    # continue
                gz += sgn * (
                    xi * np.log(yj + r)
                    + yj * np.log(xi + r)
                    - zk * np.arctan((xi * yj) / (zk * r))
                )

    return G * rho * gz * 1e5