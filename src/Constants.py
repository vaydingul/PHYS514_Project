from scipy.constants import physical_constants, m_e, m_u, hbar, c, pi, gravitational_constant

AER = 6371e3 # Average Earth radii in meters
SOLAR_MASS = 1.98847e30 # Mass of the Sun
EARTH_MASS = 5.9722e24 # Mass of the Earth
G = gravitational_constant # Gravitational constant
C = ((m_e**4) * (c**5)) / (24 * (pi ** 2) * (hbar ** 3))
D = (m_u * (m_e ** 3) * (c ** 3) * 2) / (3 * (pi ** 2) * (hbar ** 3))