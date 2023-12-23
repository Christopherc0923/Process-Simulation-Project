"""
Christopher Chan

Code for flash1 subroutine. Uses Rachford-Rice eqn.
"""

import numpy as np              # Import the module numpy and name it np
import scipy as sci             # Import the module scipy and name it sci
from scipy import optimize      # Import optimize from scipy
from Psat1 import Psat1         # Import Psat fit for Propylene
from Psat2 import Psat2         # Import Psat fit for Propane
from Psat3 import Psat3         # Import Psat fit for VCM
from Psat4 import Psat4         # Import Psat fit for Isobutane
from numpy import linalg as lin # Import linalg from scipy and name it lin

def flash1(zF, T, P):                                                                               # Subroutine called flash3, which takes in zF, T, and y1
    if 432 >= T >= 170:                                                                             # Checks if T is between 170K and 432K
        z = zF/(zF[0] + zF[1] + zF[2] + zF[3])                                                      # Creates a matrix of z 
        if (zF[0] >= 0 and zF[1] >= 0 and zF[2] >= 0 and zF[3] >= 0 and zF[0] + zF[1] + zF[2] + zF[3] > 0):  # Checks if zF components is greater than or equal to 0 and the input flow is not 0
            if P > 0:                                                                               # Checks if P is positive
                B = lambda P: (z[0]*Psat1(T)/P + z[1]*Psat2(T)/P + z[2]*Psat3(T)/P + z[3]*Psat4(T)/P) - 1   # Formula for the Bubble point pressure at T
                D = lambda P: (z[0]/Psat1(T)*P + z[1]/Psat2(T)*P + z[2]/Psat3(T)*P + z[3]/Psat4(T)*P) - 1   # Formula for the Dew point pressure at T
                Bub = sci.optimize.newton(B, 0.1)                                                   # Finds the Bubble point pressure at T
                Dew = sci.optimize.newton(D, 0.1)                                                   # Finds the Dew point pressure at T
                # print(Bub, Dew)
                if (Bub > P > Dew):     # Checks if input P is between bubble and dew pressure
                    K = [Psat1(T)/P, Psat2(T)/P, Psat3(T)/P, Psat4(T)/P]                                        # Creates an array for PsatX(T)/P where X = 1, 2, 3, and 4
                    X = [0, 0, 0, 0]    # Initalize a 4x1 matrix of 0s called X
                    y = [0, 0, 0, 0]    # Initalize a 4x1 matrix of 0s called y
                    f = lambda x: z[0] * (K[0] - 1) * (1 + x*(K[1] - 1)) * (1 + x*(K[2] - 1)) * (1 + x*(K[3] - 1)) +\
                                  z[1] * (K[1] - 1) * (1 + x*(K[0] - 1)) * (1 + x*(K[2] - 1)) * (1 + x*(K[3] - 1)) +\
                                  z[2] * (K[2] - 1) * (1 + x*(K[0] - 1)) * (1 + x*(K[1] - 1)) * (1 + x*(K[3] - 1)) +\
                                  z[3] * (K[3] - 1) * (1 + x*(K[0] - 1)) * (1 + x*(K[1] - 1)) * (1 + x*(K[2] - 1))
                                  # Rachford-Rice equation for 4 species
                    psi = sci.optimize.newton(f, 0) # Finds the root of function f, call it psi
                    V = psi * (zF[0] + zF[1] + zF[2] + zF[3])           # V/F = psi
                    L = (zF[0] + zF[1] + zF[2] + zF[3]) - V             # L = F - V
                    for i in range(4):                                  # Loops to 4
                        X[i] = z[i] / (1 + psi * (K[i] - 1))            # Rachford-Rice equation reexpressed for X[i]
                        y[i] = z[i] * K[i] / (1 + psi * (K[i] - 1))     # Rachford-Rice equation reexpressed for y[i]
                    
                    list = [X, y, V, L]                                 # Creates an array that contians liquid mole frac, vapor mole frac, V, and L. Call it list
                    return list                                         # Return list
                elif (P > Bub):
                    print("Input pressure is not in vapor liquid equilibrium, input pressure is above bubble point pressure")   # Print a string saying that P is not in VLE and return nan
                    return np.nan
                else:
                    print("Input pressure is not in vapor liquid equilibrium, input pressure is below dew point pressure")      # Print a string saying that P is not in VLE and return nan
                    return np.nan
            else:
                print("Input pressure is not positive"),                               # Print a string saying that P is not valid and return nan
                return np.nan

        else:
            print("The input feed component molar flow rate is not valid.")         # Print a string saying that zF is not in VLE and return nan
            return np.nan

    else:
        print("The input temperature is not within the range.")                     # Print a string saying that T is not in VLE and return nan
        return np.nan


