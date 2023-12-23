"""
Christopher Chan

Code for IB fit using lagrange polynomials
"""
# Christopher Chan Hours Spent: 20
import numpy as np  # Import the module numpy and name it np so we can use the exp function

# Two lagrange polynomials were constructed for better approximation of the vapor pressure, these lagrange polynomials are piecewise at T = 300K
def Psat4(T):               # Creates a subroutine Psat4 which takes in x and returns the vapor pressure at T (if possible)
    if  300 >= T >= 170:    # Checks if the input is between 170 to 300 for the first lagrange polynomial 
        x = 1 / T           # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 170 to T = 300
        Poly = -12944430449216896.0*x**6 + 341351153523758.30104118104118104*x**5 - 3736321842954.1413686034130478575*x**4 + 21737180182.86316536763647874759*x**3 - 71023633.255350749495589001761841*x**2 + 121165.12503156485334510025868051*x - 79.946337506279189489066032275909
        # The Lagrange polynomials looks relatively linear for ln(P*) vs T, which means that here is a linear relationship between the two.
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    if 432 >= T > 300:      # Checks if the input is between 300 to 420 for the second lagrange polynomial 
        x = 1/T             # Converts the input temperature to 1/T
        
        # Lagrange polynomial evaluates ln(P*) from T = 300.01 to T = 432
        Poly = -16947712453680.0*x**5 + 346946919810.62518403852769177847*x**4 - 2719539756.152731260749914000688*x**3 + 10316922.760192036965370943699117*x**2 - 21642.847697927653078775369797042*x + 23.68514267943416036005045292971
        VaporPressure = np.exp(Poly)    # The value evaluated from the lagrange polynomial is exponentiated with natural base to calculate the vapor pressure
        return VaporPressure            # Returns the vapor pressure in bars

    else:                               # Input temperature is not in the range, returns nan
        return np.nan 
