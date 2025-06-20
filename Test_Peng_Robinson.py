import numpy as np
from scipy.optimize import fsolve
import math

def equation(Z,A,B,):
    return Z**3 - (1 - B) * Z**2 + (A -2*B-3*B**2)*Z - (A*B - B**2 - B**3) 

#DATA
comp = ["CH4" , "C2H6", "H2"]
T_c_K = [190.6, 305.3, 33.2]
P_c_bar = [45.99, 48.72, 12.8]
omega = [0.011, 0.099, -0.218] #acentric factor
x = [0.4,0.4, 0.2] #molar fraction
k_matrix = np.array([[0, 0.01, 0.02],
                     [0.01, 0, 0.015],
                     [0.02, 0.015, 0]])

initial_guess = 1.0 #INITIAL GUESS for Z

#constants
R = 0.08314

#input
T = 300 #K
P = 50 #bar

print("Welcome! The current parameters for this simulation are:")
print("Component of the mixture and their molar composition:")
for i in range(len(comp)):
    print(f"{comp[i]}\t{x[i]}")
print("The parameters of the simulation are:")
print("T = ", T, " K")
print("P = ",P, "bar")
print("------------------------------------------------------------------")
#initialization PR parameters
a = [0.0] * len(comp)
b = [0.0] * len(comp)

for i in range(len(comp)):
    a[i] = (0.45724 * (R*T_c_K[i])**2 / P_c_bar[i] )*( 1 + (0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i]**2)*(1-math.sqrt(T/T_c_K[i])))**2
    b[i] = 0.07780 * (R*T_c_K[i])/ P_c_bar[i]  
print("Coefficents for the calculation:")
print("a = ",a)
print("b = ",b)

a_m = 0.0
for i in range(len(comp)):
    for j in range(len(comp)):
        k = k_matrix[i,j]
        a_m = a_m + (x[i] * x[j] * math.sqrt(a[i] * a[j]) * (1 - k))
        #print("Indexes: ", i,j)
print("a_m = " ,a_m)

b_m = 0.0
for i in range(len(comp)):
    b_m = b_m + x[i] * b[i]

print("b_m = " ,b_m)

A = (a_m * P)/ (R**2 * T**2)
B = (b_m * P)/(R*T)
print("A = ", A)
print("B = ", B)
print("-------------------------------------------------------------------")
print("Results:")
#compressibility calculation
Z = fsolve(equation,initial_guess,args = (A,B))
Z = Z[0]
print("Z = ", round(Z,4))

#fugacity coeff calculation
fug_coef = [0.0] * len(comp)

for i in range(len(comp)):
    sumj = 0.0
    for j in range(len(comp)):
        k = k_matrix[i,j]
        sumj = sumj + (x[j] * math.sqrt(a[i] * a[j])*(1-k))
        #print(i,j)    
    fug_coef[i] = math.exp((b[i]/b_m)*(Z-1)- math.log(Z - B)- (A/(2*math.sqrt(2)*B)) *((b[i]/b_m) - (2/a_m) * sumj )* math.log((Z+(1+math.sqrt(2)) * B)/(Z+(1 - math.sqrt(2)) * B)))
print("fugacity coefficients = ", fug_coef)

#fugacity calculation
fug = [0.0] * len(comp)        
for i in range(len(comp)):
    fug[i] = fug_coef[i] * x[i] * P
print("Fugacity of the component of the mixture (bar):")
for i in range(len(comp)):
    print(f"{comp[i]}\t{round(fug[i],3)}")
