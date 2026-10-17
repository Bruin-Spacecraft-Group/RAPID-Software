"""
Determines the attitude of the RAPID-0 satellite with two vectors from sensors.
"""

import ulab.numpy as np
from quaternion import Quaternion

# Status codes for each case
SUCCESS = 0  # Success: attitude estimated and returned successfully
ANTI_PARALLEL = 1  # Anti-parallel vectors: estimated attitude with 180 degree rotation
COLLINEAR = 2  # Collinear vectors: failure, impossible to estimate attitude
SINGULAR = 3  # Singular: failure from insufficient information to estimate attitude
NORM_ERR = 4  # Normalization error: failure from prevented division by zero

# pylint: disable=too-many-locals
def triad_algorithm(
    r1: np.ndarray, r2: np.ndarray, b1: np.ndarray, b2: np.ndarray
) -> tuple[Quaternion, int]:
    """
    Calculates the attitude using an optimized TRIAD algorithm based off of:
    'Fast Quaternion Attitude Estimation From Two Vector Measurements' (Markley, 2002).

    Args:
        r1 (np.ndarray): The first (more accurate) reference vector in inertial frame 
                         (e.g. Sun vector from environment model).
        r2 (np.ndarray): The second (less accurate) reference vector in inertial frame 
                         (e.g. Magnetic field from environment model).
        b1 (np.ndarray): The measurement of the first vector in the body frame 
                         (e.g. Sun sensor reading).
        b2 (np.ndarray): The measurement of the second vector in the body frame 
                         (e.g. Magnetometer reading).

    Returns:
        tuple[Quaternion, int]:
            Quaternion: The estimated attitude Quaternion from the body frame to the inertial frame
                        (identity Quaternion if calculation is impossible with given input vectors)
            int: A status code indicating the result.
    """

    # Normalize input vectors
    r1_unit = r1 / np.linalg.norm(r1)
    r2_unit = r2 / np.linalg.norm(r2)
    b1_unit = b1 / np.linalg.norm(b1)
    b2_unit = b2 / np.linalg.norm(b2)

    # Calculate intermediate cross products (not normalised)
    r3 = np.cross(r1_unit, r2_unit)
    b3 = np.cross(b1_unit, b2_unit)

    # Check for collinearity (which gives a zero cross product, failure case)
    if np.linalg.norm(b3) < 1e-9 or np.linalg.norm(r3) < 1e-9:
        return Quaternion(), COLLINEAR

    # Check if b1 and r1 are anti-parallel (to prevent division by zero)
    dot_b1_r1 = np.dot(b1_unit, r1_unit)
    one_plus_dot = 1.0 + dot_b1_r1

    # If b1 and r1 are anti-parallel, return 180 degree rotation quaternion from a vector perp to r1
    if one_plus_dot < 1e-9:
        c = np.cross(r2_unit, b2_unit)  # Construct a vector that maps r2 to b2
        # Project c onto the plane orthogonal to r1_unit
        c_proj = c - np.dot(c, r1_unit) * r1_unit
        c_proj_norm = np.linalg.norm(c_proj)

        # Check if c is parallel to r1, and pick a perpendicular
        if c_proj_norm < 1e-9:
            # Find the index of the smallest absolute component of r1 to maximize orthogonality
            min_axis = np.argmin(np.abs(r1_unit))
            # Construct the corresponding basis vector
            basis = np.zeros(3)
            basis[min_axis] = 1.0
            axis = np.cross(r1_unit, basis)  # Construct a vector perpendicular to r1
            axis /= np.linalg.norm(axis)
        else:
            axis = c_proj / c_proj_norm

        return Quaternion(0.0, *axis), ANTI_PARALLEL

    # Eq. (31) from the Markley 2002 paper
    mu = one_plus_dot * np.dot(b3, r3) - np.dot(b1_unit, r3) * np.dot(r1_unit, b3)

    b1_plus_r1 = b1_unit + r1_unit
    nu = np.dot(b1_plus_r1, np.cross(b3, r3))  # Eq. (32) from the Markley paper

    rho = np.sqrt(mu**2 + nu**2)  # Eq. (33) from the Markley paper

    # If rho is near zero, then mu and nu are also near zero
    if rho < 1e-9:
        return Quaternion(), SINGULAR

    if mu >= 0:  # Eq. (35a) from the Markley paper
        rho_plus_mu = rho + mu

        # Calculate the square root term in the denominator of the scalar multiple
        sqrt_term = rho * rho_plus_mu * one_plus_dot
        if sqrt_term < 1e-9:  # Prevent division by zero
            return Quaternion(), NORM_ERR

        # The scalar multiple
        mult = 0.5 / np.sqrt(sqrt_term)

        # The scalar part of the quaternion
        w = rho_plus_mu * one_plus_dot

        # The vector part of the quaternion
        v = rho_plus_mu * np.cross(b1_unit, r1_unit) + nu * b1_plus_r1

        # Instantiate the estimated attitude Quaternion
        q = Quaternion(w * mult, v[0] * mult, v[1] * mult, v[2] * mult)

    else:  # Eq. (35b) from the Markley paper
        rho_minus_mu = rho - mu

        # Calculate the square root term in the denominator of the scalar multiple
        sqrt_term = rho * rho_minus_mu * one_plus_dot
        if sqrt_term < 1e-9:  # Prevent division by zero
            return Quaternion(), NORM_ERR

        # The scalar multiple
        mult = 0.5 / np.sqrt(sqrt_term)

        # The scalar part of the quaternion
        w = nu * one_plus_dot

        # The vector part of the quaternion
        v = nu * np.cross(b1_unit, r1_unit) + rho_minus_mu * b1_plus_r1

        # Instantiate the estimated attitude Quaternion
        q = Quaternion(w * mult, v[0] * mult, v[1] * mult, v[2] * mult)

    return q, SUCCESS
