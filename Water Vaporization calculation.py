import math
from scipy.optimize import fsolve
import numpy as np

#inizialization for controls
mode = "A" 
T_C = -300.0
P_in = -1 

#DIPPR Coefficients
A = 73.649
B = -7258.2
C = -7.3037
D = 4.17E-6
E = 2

def P_vap():
    P_Pa = math.exp(A + B/T + C * math.log(T)+D * T**E)
    return P_Pa

# Define the function to solve
def equation(T, A, B, C, D, E, P):
    return np.exp(A + B/T + C*np.log(T) + D*T**E) - P

def T_vap(P):
    initial_guess = 350.0
    # Solve for T using fsolve
    T_vap = fsolve(equation, initial_guess, args=(A, B, C, D, E, P))
    return T_vap

#P[Pa] = exp(A + B/T + ClnT + DT^E), T in [K]
#T [K] = (a+blnP)/(1+clnP+d(lnP)^2), P in [Pa]

def Antoine():
    A_A = 6.20963
    B_A = 2354.731
    C_A = 7.559
    P_A = 10**(A_A-(B_A/(T+C_A))) *1E5
    return P_A

while mode != "T" and mode != "P" :
    mode = input("Do you want to calculate Pressure of vaporization or Temperature for a given Pressure? (T or P):" )

if mode == "P":
    print("Entering in the pressure calculation mode.")
    while T_C < -273.15: 
        T_C = float(input("Enter the Temperature (°C): "))
    T = T_C + 273.15
    P_Pa = P_vap()
    print("The vaporization pressure for this Temperature is: ", round(P_Pa,2), " Pa")
    if 293 < T < 343:
        P_A = Antoine()
        print("Since the T is between 20 and 70 °C, you may consider also the Pressure given by the Antoine Equation: ", round(P_A,2)," Pa")
else:
    print("Entering in the temperature calculation mode.")
    while P_in < 0: 
        P_in = float(input("Enter the Pressure (Pa): "))
    T_vap = T_vap(P_in)
    T_vap_C = T_vap[0] - 273.15
    print("The numerical approximate vaporization Temperature for this Pressure is: ", round(T_vap_C,2), " °C")
