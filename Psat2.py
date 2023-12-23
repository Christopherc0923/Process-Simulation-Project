"""
Christopher Chan

Code for P fit using lagrange polynomials
"""
# Christopher Chan Hours Spent: 20
import numpy as np  # Import the module numpy and name it np so we can use the exp function

# Two lagrange polynomials were constructed for better approximation of the vapor pressure, these lagrange polynomials are piecewise at T = 300K
def Psat2(T):               # Creates a subroutine Psat2 which takes in x and returns the vapor pressure at T (if possible)
    if  432 >= T >= 170:    # Checks if the input is between 170 to 432 
        x = 1 / T           # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 170 to T = 432
        Poly = -9.55703314e+20*x**8 + 3.21181073e+19*x**7 + -4.67666254e+17*x**6 + 3.85043573e+15*x**5 + -1.95817511e+13*x**4 + 6.28685199e+10*x**3 + -1.24135231e+08*x**2 + 1.35164594e+05*x + -5.52481146e+01
        # The Lagrange polynomials looks relatively linear for ln(P*) vs T, which means that here is a linear relationship between the two.
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    else:                               # Input temperature is not in the range, returns nan
        return np.nan 
    
