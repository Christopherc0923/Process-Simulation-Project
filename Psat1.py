"""
Christopher Chan

Code for PM fit using lagrange polynomials
"""
# Christopher Chan Hours Spent: 20
import numpy as np  # Import the module numpy and name it np so we can use the exp function

# Two lagrange polynomials were constructed for better approximation of the vapor pressure, these lagrange polynomials are piecewise at T = 300K
def Psat1(T):               # Creates a subroutine Psat1 which takes in x and returns the vapor pressure at T (if possible)
    if  432 >= T >= 170:    # Checks if the input is between 170 to 432 
        x = 1 / T           # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 170 to T = 432
        # Poly = -9.66892959e+17*x**7 + 3.08095861e+16*x**6 -4.20210369e+14*x**5 + 3.18406937e+12*x**4 + -1.44982687e+10*x**3 + 3.96709034e+07*x**2 + -6.25327418e+04*x + 4.91128752e+01
        Poly = (9.23338745e+19*x**7 + -2.92235493e+18*x**6 + 3.94456853e+16*x**5 + -2.94373231e+14*x**4 + 1.31185315e+12*x**3 + -3.49139546e+09*x**2 + 5.13658189e+06*x + -3.21726058e+03)
        # The Lagrange polynomials looks relatively linear for ln(P*) vs T, which means that here is a linear relationship between the two.
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    else:                               # Input temperature is not in the range, returns nan
        return np.nan 
