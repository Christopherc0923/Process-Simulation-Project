"""
Christopher Chan

Code for flash2 subroutine. Use steepest descent with simple line search.
"""
import numpy as np              # Import the module numpy and name it np
import scipy as sci             # Import the module scipy and name it sci
from scipy import optimize      # Import optimize from scipy
from Psat1 import Psat1         # Import Psat fit for Propylene
from Psat2 import Psat2         # Import Psat fit for Propane
from Psat3 import Psat3         # Import Psat fit for VCM
from Psat4 import Psat4         # Import Psat fit for Isobutane
from numpy import linalg as lin # Import linalg from scipy and name it lin

def flash2(zF, T, L):                                                                               # Subroutine called flash3, which takes in zF, T, and y1
    if 432 >= T >= 170:                                                                             # Checks if T is between 170K and 432K
        if (zF[0] > 0 and zF[1] > 0 and zF[2] > 0 and zF[3] > 0):                                   # Checks if zF components is positive
            z = zF/(zF[0] + zF[1] + zF[2] + zF[3])                                                  # Creates a matrix of z
            F = zF[0] + zF[1] + zF[2] + zF[3]                                                       # Evaluates the feed flow rate
            B = lambda P: (zF[0]/F*Psat1(T)/P + zF[1]/F*Psat2(T)/P + zF[2]/F*Psat3(T)/P + zF[3]/F*Psat4(T)/P) - 1 # Formula for the Bubble point pressure at T
            D = lambda P: (zF[0]/F/Psat1(T)*P + zF[1]/F/Psat2(T)*P + zF[2]/F/Psat3(T)*P + zF[3]/F/Psat4(T)*P) - 1 # Formula for the Dew point pressure at T
            Bub = sci.optimize.newton(B, 0.1)                                                   # Finds the Bubble point pressure at T
            Dew = sci.optimize.newton(D, 0.1)                                                   # Finds the Dew point pressure at T
            if F > L and L > 0:                                                                     # Checks if liquid flow rate less than feed flow rate and if liquid flow rate is positive
                V = F - L                           # Evaluates V by subtracting F by L
                x = [.25, .25, .25, .25, (Bub + Dew) / 2]    # Creates a 5x1 matrix of the inital guess for x1, x2, x3, x4, and P. Call it x
                TOL = 1e-8                                  # Sets the Tolerance to 1e-8
                j = 0                                       # Set step size exponent to 0
                for counter in range(20000):                # Loops to 20000
                    f1 = - Psat1(T)*x[0] - Psat2(T)*x[1] - Psat3(T)*x[2] - Psat4(T)*x[3] + x[4]     # Defines and evaluates function 1: P - Psat1(T)*x1 - Psat2(T)*x2 - Psat3(T)*x3 - Psat4(T)*x4 
                    f2 = z[0]*F*x[4] - Psat1(T)*x[0]*V - x[4]*x[0]*L                                # Defines and evaluates function 2: zF1*P - Psat1(T)*x1*V - x1*L*P
                    f3 = z[1]*F*x[4] - Psat2(T)*x[1]*V - x[4]*x[1]*L                                # Defines and evaluates function 3: zF2*P - Psat2(T)*x2*V - x2*L*P
                    f4 = z[2]*F*x[4] - Psat3(T)*x[2]*V - x[4]*x[2]*L                                # Defines and evaluates function 4: zF3*P - Psat3(T)*x3*V - x3*L*P
                    f5 = z[3]*F*x[4] - Psat4(T)*x[3]*V - x[4]*x[3]*L                                # Defines and evaluates function 5: zF4*P - Psat4(T)*x4*V - x4*L*P
                    def g(x):                                                                       # Subroutine called g, which takes in x and returns f1^2 + f2^2 + f3^2 + f4^2 + f5^2
                        f1 = x[4] - Psat1(T)*x[0] - Psat2(T)*x[1] - Psat3(T)*x[2] - Psat4(T)*x[3]
                        f2 = z[0]*F*x[4] - Psat1(T)*x[0]*V - x[4]*x[0]*L
                        f3 = z[1]*F*x[4] - Psat2(T)*x[1]*V - x[4]*x[1]*L
                        f4 = z[2]*F*x[4] - Psat3(T)*x[2]*V - x[4]*x[2]*L
                        f5 = z[3]*F*x[4] - Psat4(T)*x[3]*V - x[4]*x[3]*L
                        return f1**2 + f2**2 + f3**2 + f4**2 + f5**2
                    
                    def gradf(x):   # Subroutine called gradf, which takes in x and returns the gradient of f1^2 + f2^2 + f3^2 + f4^2 + f5^2
                        return 2.0*f1*np.array([-Psat1(T), -Psat2(T), -Psat3(T), -Psat4(T), 1]) +\
                                2.0*f2*np.array([-Psat1(T)*V - L*x[4], 0, 0, 0, z[0]*F-L*x[0]]) +\
                                2.0*f3*np.array([0, -Psat2(T)*V - L*x[4], 0, 0, z[1]*F-L*x[1]]) +\
                                2.0*f4*np.array([0, 0, -Psat3(T)*V - L*x[4], 0, z[2]*F-L*x[2]]) +\
                                2.0*f5*np.array([0, 0, 0, -Psat4(T)*V - L*x[4], z[3]*F-L*x[3]])
                    d = -gradf(x)           # Search direction is the negative of the gradient
                    t = 0.5**j              # Step size
                    if g(x + t*d) > g(x):   # Checks if g(x + t*d) is greater than g(x)
                        j = j + 1           # If true, increase j by 1, which will decrease the step size
                    
                    elif t * lin.norm(d)/lin.norm(x) < TOL: # Checks if stop condition has been reached
                        X = [x[0], x[1], x[2], x[3]]        # Compiles an array of the liquid mole fraction
                        y = [x[0]*Psat1(T)/x[4], x[1]*Psat2(T)/x[4], x[2]*Psat3(T)/x[4], x[4]*Psat4(T)/x[4]]    # Compiles an array of the vapor mole fraction
                        P = x[4]                                                                                # Sets P as x[4]
                        list = [X, y, V, P]                                                                     # Compiles an array, called list, that contians the liquid mole fracs, vapor mole fracs, pressure, vapor molar flow rate, and pressure
                        return list                                                                             # Returns the list
                        break                                                                                   # Stops the loop
                    else:                                                                                       # Sets new guess of x to x + t*d
                        x = x + t*d
                print ("Max iterations Reached")                                                             # Returns max iteration is reached 
                X = [x[0], x[1], x[2], x[3]]        # Compiles an array of the liquid mole fraction
                y = [x[0]*Psat1(T)/x[4], x[1]*Psat2(T)/x[4], x[2]*Psat3(T)/x[4], x[4]*Psat4(T)/x[4]]    # Compiles an array of the vapor mole fraction
                P = x[4]                                                                                # Sets P as x[4]
                list = [X, y, V, P]                                                                     # Compiles an array, called list, that contians the liquid mole fracs, vapor mole fracs, pressure, vapor molar flow rate, and pressure
                return list 
               
            else:
                print("The liquid molar flow rate is not valid")                                    # Prints a string saying that V is not valid and return nan
                return np.nan
        else:
            print("The input feed component molar flow rate is not valid.")                         # Prints a string saying that zF is not valid and return nan
            return np.nan
    else:
        print("The input temperature is not within the VLE.")                                       # Prints a string saying that T is not valid and return nan
        return np.nan

# print(flash2(np.array([0.25, 0.25, 0.25, 0.25]), 310, 0.5031898890192623))