"""
Christopher Chan

Code for CCflash subroutine
"""
import numpy as np              # Import the module numpy and name it np
import scipy as sci             # Import the module scipy and name it sci
from scipy import optimize      # Import optimize from scipy
from Psat1 import Psat1         # Import Psat fit for Propylene
from Psat2 import Psat2         # Import Psat fit for Propane
from Psat3 import Psat3         # Import Psat fit for VCM
from Psat4 import Psat4         # Import Psat fit for Isobutane
from flash1 import flash1       # Import flash1 subroutine
from flash2 import flash2       # Import flash2 subroutine
from flash3 import flash3       # Import flash3 subroutine


def CCflash(zF, T, P, L, y1):   # Subroutine called CCflash, which takes in zF, T, P, L, and y1
    if 432 >= T >= 170:                                                                                         # Checks if T is between 170K and 432K
        F = zF[0] + zF[1] + zF[2] + zF[3]                                                                       # Sums the feed compenent molar flow 
        if (zF[0] >= 0 and zF[1] >= 0 and zF[2] >= 0 and zF[3] >= 0 and zF[0] + zF[1] + zF[2] + zF[3] > 0):     # Checks if zF components is greater than or equal to 0 and the input flow is not 0                          
            z = zF/(zF[0] + zF[1] + zF[2] + zF[3])
            B = lambda P: (z[0]*Psat1(T)/P + z[1]*Psat2(T)/P + z[2]*Psat3(T)/P + z[3]*Psat4(T)/P) - 1           # Formula for the Bubble point pressure at T
            D = lambda P: (z[0]/Psat1(T)*P + z[1]/Psat2(T)*P + z[2]/Psat3(T)*P + z[3]/Psat4(T)*P) - 1           # Formula for the Dew point pressure at T
            Bub = sci.optimize.newton(B, 0.1)                                                                   # Finds the Bubble point pressure at T
            Dew = sci.optimize.newton(D, 0.1)                                                                   # Finds the Dew point pressure at T
            if P > 0 and Bub > P > Dew:                                                                         # Checks if P is positive and if the pressure is in vapor-liquid equilibrium
                print("zF, T, and P are valid, use flash1")
                list = [flash1(zF, T, P)[0], flash1(zF, T, P)[1], flash1(zF, T, P)[2], flash1(zF, T, P)[3], P]  # Creates an array that contians liquid mole frac, vapor mole frac, V, L, and P. Call it list
                return list                                                                                     # Returns list
            elif(F > L and L > 0):
                print("P not positive or not in VLE, so flash2 was used")       
                list = [flash2(zF, T, L)[0], flash2(zF, T, L)[1], flash2(zF, T, L)[2], L, flash2(zF, T, L)[3]]  # Creates an array that contians liquid mole frac, vapor mole frac, V, L, and P. Call it list
                return list                                                                                     # Returns list
            elif(1 >= y1 >= 0):
                print("P not positive or not in VLE and L was not valid, so flash3 was used")
                list = flash3(zF, T, y1)
                return list
            else:
                print("No method work")
                return np.nan
        else:
            print("The input feed component molar flow rate is not valid.")         # Print a string saying that zF is not in VLE and return nan
            return np.nan
    else:
        print("The input temperature is not within the range.")                     # Print a string saying that T is not in VLE and return nan
        return np.nan
