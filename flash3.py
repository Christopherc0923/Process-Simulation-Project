"""
Christopher Chan

Code for flash3 subroutine
"""
import numpy as np              # Import the module numpy and name it np
import scipy as sci             # Import the module scipy and name it sci
from scipy import optimize      # Import optimize from scipy
from Psat1 import Psat1         # Import Psat fit for Propylene
from Psat2 import Psat2         # Import Psat fit for Propane
from Psat3 import Psat3         # Import Psat fit for VCM
from Psat4 import Psat4         # Import Psat fit for Isobutane
from numpy import linalg as lin # Import linalg from scipy and name it lin


def flash3(zF, T, y1):                                                                              # Subroutine called flash3, which takes in zF, T, and y1
    if 432 >= T >= 170:                                                                             # Checks if T is between 170K and 432K
        if (zF[0] > 0 and zF[1] > 0 and zF[2] > 0 and zF[3] > 0):                                   # Checks if zF components is positive
            if 1 >= y1 >= 0:                                                                        # Checks if y1 is between 0 and 1
                
                b = np.array([[0], [0], [0], [0], [0], [0], [0], [0], [0], [0]])                    # Initial search direction as a 10x1 matrix of 0 called b
                F = zF[0] + zF[1] + zF[2] + zF[3]                                                   # Evaluate the value of F
                B = lambda P: (zF[0]/F*Psat1(T)/P + zF[1]/F*Psat2(T)/P + zF[2]/F*Psat3(T)/P + zF[3]/F*Psat4(T)/P) - 1 # Formula for the Bubble point pressure at T
                D = lambda P: (zF[0]/F/Psat1(T)*P + zF[1]/F/Psat2(T)*P + zF[2]/F/Psat3(T)*P + zF[3]/F/Psat4(T)*P) - 1 # Formula for the Dew point pressure at T
                Bub = sci.optimize.newton(B, 0.1)                                                   # Finds the Bubble point pressure at T
                Dew = sci.optimize.newton(D, 0.1)                                                   # Finds the Dew point pressure at T
                xk = np.array([[.25], [.25], [.25], [.25], [.25], [.25], [.25], [(Bub + Dew)/ 2], [0.5*F], [0.5*F]]) # Initial guess for x1, x2, x3, x4, y2, y3, y4, P, L, V as a 10x1 matrix called xk 
                Tol = 1e-8                                                                          # Sets the Tolerance to 1e-8
                for counter in range(10000):                                                        # Loops to 10000
            
                    f1 = 1 - xk[0, 0] - xk[1, 0] - xk[2, 0] - xk[3, 0]                              # Defines and evaluates function 1: 1 - x1 - x2 - x3 - x4
                    f2 = 1 - y1 - xk[4, 0] - xk[5, 0] - xk[6, 0]                                    # Defines and evaluates function 2: 1 - y1 - y2 - y3 - y4
                    f3 = zF[0] - y1*xk[9, 0] - xk[0, 0]*xk[8, 0]                                    # Defines and evaluates function 3: zF1 - y1*V - x1*L
                    f4 = zF[1] - xk[4, 0]*xk[9, 0] - xk[1, 0]*xk[8, 0]                              # Defines and evaluates function 4: zF2 - y2*V - x2*L
                    f5 = zF[2] - xk[5, 0]*xk[9, 0] - xk[2, 0]*xk[8, 0]                              # Defines and evaluates function 5: zF3 - y3*V - x3*L
                    f6 = zF[3] - xk[6, 0]*xk[9, 0] - xk[3, 0]*xk[8, 0]                              # Defines and evaluates function 6: zF4 - y4*V - x4*L
                    f7 = xk[0, 0] - y1*xk[7, 0]/Psat1(T)                                            # Defines and evaluates function 7: x1 - y1*P/Psat1(T), Raoult's Law
                    f8 = xk[1, 0] - xk[4, 0]*xk[7, 0]/Psat2(T)                                      # Defines and evaluates function 8: x2 - y2*P/Psat2(T), Raoult's Law
                    f9 = xk[2, 0] - xk[5, 0]*xk[7, 0]/Psat3(T)                                      # Defines and evaluates function 9: x3 - y3*P/Psat3(T), Raoult's Law
                    f10 = xk[3, 0] - xk[6, 0]*xk[7, 0]/Psat4(T)                                     # Defines and evaluates function 10: x4 - y4*P/Psat4(T), Raoult's Law
                    f = np.array([[f1],[f2],[f3],[f4],[f5],[f6],[f7],[f8],[f9],[f10]])              # Creates a 10x1 matrix of the functions called f
            
                    J = np.array([[-1, -1, -1, -1, 0, 0, 0, 0, 0, 0],                               # Define the Jacobian using function f1 - f10. Takes the parital derivative of each function with respect to either x1, x2, x3, x4, y2, y3, y4, P, L, or V
                                  [0, 0, 0, 0, -1, -1, -1, 0, 0, 0],
                                  [-xk[8, 0], 0, 0, 0, 0, 0, 0, 0, -xk[0, 0], -y1],
                                  [0, -xk[8, 0], 0, 0, -xk[9, 0], 0, 0, 0, -xk[1, 0], -xk[4, 0]],
                                  [0, 0, -xk[8, 0], 0, 0, -xk[9, 0], 0, 0, -xk[2, 0], -xk[5, 0]],
                                  [0, 0, 0, -xk[8, 0], 0, 0, -xk[9, 0], 0, -xk[3, 0], -xk[6, 0]],
                                  [1, 0, 0, 0, 0, 0, 0, -y1/Psat1(T), 0, 0],
                                  [0, 1, 0, 0, -xk[7, 0]/Psat2(T), 0, 0, -xk[4, 0]/Psat2(T), 0, 0],
                                  [0, 0, 1, 0, 0, -xk[7, 0]/Psat3(T), 0, -xk[5, 0]/Psat3(T), 0, 0],
                                  [0, 0, 0, 1, 0, 0, -xk[7, 0]/Psat4(T), -xk[6, 0]/Psat4(T), 0, 0]])
                                  
                    b = np.linalg.solve(J, f)                                                       # Solves the matrix called b which satisfies J@b = f
                    if lin.norm(b, ord = np.inf) < Tol:                                             # Stop condition that checks if the infinite norm of b is less than Tol
                        x = np.array([xk[0, 0], xk[1, 0], xk[2, 0], xk[3, 0]])                      # Compiles an array of the liquid mole fraction
                        y = np.array([y1, xk[4, 0], xk[5, 0], xk[6, 0]])                            # Compiles an array of the vapor mole fraction
                        if xk[9, 0] > 0 and xk[8, 0] > 0 and xk[7, 0]:
                            V = xk[9, 0]                                                                # Defines V from the xk matrix 
                            L = xk[8, 0]                                                                # Defines L from the xk matrix
                            P = xk[7, 0]                                                                # Defines P from the xk matrix
                            list = [x, y, V, L, P]                                                      # Compiles an array, called list, that contians the liquid mole fracs, vapor mole fracs, pressure, liquid molar flow rate, and vapor molar flow rate
                            return list                                                                 # Returns the list
                            break                                                                       # Stops the loop
                        else:
                            print("flash3 did not work with the inputs")
                            return np.nan
                            break
                    else:
                        xk = xk - b                                                                 # If the stop conditions has not been satisfied, set the new xk to xk - b
                return("Max iterations Reached", xk)                                                # Returns max iteration is reached and the current guess
            else:
                print("The input vapor fraction is not between 0 and 1.")                           # Return a string saying that y1 is not valid and return nan
                return np.nan
        else:
            print("The input feed component molar flow rate is not valid.")                         # Return a string saying that zF is not valid and return nan
            return np.nan
    else:
        print("The input temperature is not the desired range.")                                    # Return a string saying that T is not valid and return nan
        return np.nan

