"""
Unit Conversion Functions for Nuclear Fission Plants Project
Utility module for common unit conversions
"""

from scipy.special import j0, j1


def psia_to_pa(pressure_psia):
    """
    Convert pressure from psia (pounds per square inch absolute) to Pa (Pascals).
    
    Parameters:
    - pressure_psia: float, pressure in psia
    
    Returns:
    - float, pressure in Pa
    """
    return pressure_psia * 6894.757


def lbm_per_hr_to_kg_per_s(mass_flow_lbm_hr):
    """
    Convert mass flow rate from lbm/hr (pounds mass per hour) to kg/s (kilograms per second).
    
    Parameters:
    - mass_flow_lbm_hr: float, mass flow rate in lbm/hr
    
    Returns:
    - float, mass flow rate in kg/s
    """
    return mass_flow_lbm_hr * 0.000125998


def fahrenheit_to_celsius(temp_fahrenheit):
    """
    Convert temperature from °F (Fahrenheit) to °C (Celsius).
    
    Parameters:
    - temp_fahrenheit: float, temperature in °F
    
    Returns:
    - float, temperature in °C
    """
    return (temp_fahrenheit - 32) * 5 / 9


def inches_to_meters(length_inches):
    """
    Convert length from in (inches) to m (meters).
    
    Parameters:
    - length_inches: float, length in inches
    
    Returns:
    - float, length in meters
    """
    return length_inches * 0.0254


def square_feet_to_square_meters(area_sq_ft):
    """
    Convert area from ft² (square feet) to m² (square meters).
    
    Parameters:
    - area_sq_ft: float, area in square feet
    
    Returns:
    - float, area in square meters
    """
    return area_sq_ft * 0.092903


def bessel_j0(x):
    """
    Bessel function of the first kind of order 0.
    
    Parameters:
    - x: float or array, argument of the Bessel function
    
    Returns:
    - float or array, J0(x)
    """
    return j0(x)


def bessel_j1(x):
    """
    Bessel function of the first kind of order 1.
    
    Parameters:
    - x: float or array, argument of the Bessel function
    
    Returns:
    - float or array, J1(x)
    """
    return j1(x)
