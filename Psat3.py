"""
Christopher Chan

Code for VCM fit using lagrange polynomials
"""
# Christopher Chan Hours Spent: 20
import numpy as np  # Import the module numpy and name it np so we can use the exp function

# Two lagrange polynomials were constructed for better approximation of the vapor pressure, these lagrange polynomials are piecewise at T = 300K
def Psat3(T):               # Creates a subroutine Psat3 which takes in x and returns the vapor pressure at T (if possible)
    if  300 >= T >= 170:    # Checks if the input is between 170 to 300 for the first lagrange polynomial 
        x = 1 / T           # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 170 to T = 300
        Poly = 59624381631218560.0*x**6 - 1572940844017440.4276967476967477*x**5 + 17186118123599.868132815332815333*x**4 - 99533958476.633880244693578026911*x**3 + 322059927.57915469540336207002874*x**2 - 554322.09649597389894056560723227*x + 400.93120217951523784857118190452
        # The Lagrange polynomials looks relatively linear for ln(P*) vs T, which means that here is a linear relationship between the two.
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    if 432 >= T > 300:      # Checks if the input is between 300 to 420 for the second lagrange polynomial 
        x = 1/T             # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 300.01 to T = 432
        Poly = 3583795805511901184.0*x**6 - 60669166353504240.470415253820827*x**5 + 426937662920148.11761681327501761*x**4 - 1598815694969.7700535451919559247*x**3 + 3360854921.3027941413446470206016*x**2 - 3763076.9469543930195860020421424*x + 1759.9146315940153215282214250223
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    else:                               # Input temperature is not in the range, returns nan
        return np.nan 
