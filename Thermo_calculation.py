#functions

def aprox(input_to_aprox):
    decimal_positions = 3
    output_aprox = round(input_to_aprox,decimal_positions)
    return output_aprox
    
def MM_tot():
    MM_tot = 0.0
    for i in range(len(comp)):
        MM_tot = MM_tot + comp[i]/100 * MM[i] 
    return MM_tot

def calc_Cp_mix():
    Acp_tot = 0.0
    Bcp_tot = 0.0
    Ccp_tot = 0.0
    Dcp_tot = 0.0
    Ecp_tot = 0.0
    for i in range(len(comp)):
        Acp_tot = Acp_tot + comp[i]/100 * Acp[i]
    for i in range(len(comp)):
        Bcp_tot = Bcp_tot + comp[i]/100 * Bcp[i]
    for i in range(len(comp)):
        Ccp_tot = Ccp_tot + comp[i]/100 * Ccp[i]
    for i in range(len(comp)):
        Dcp_tot = Dcp_tot + comp[i]/100 * Dcp[i]
    for i in range(len(comp)):
        Ecp_tot = Ecp_tot + comp[i]/100 * Ecp[i]
    Cp_mix = 0.0
    Cp_mix = Acp_tot + Bcp_tot * T + Ccp_tot * T**2 + Dcp_tot * T**3 + Ecp_tot * T**-2
    return Cp_mix

def calc_Cp():
    Cp = []
    for i in range(len(comp)):
        Cp.append(Acp[i] + Bcp[i] * T + Ccp[i] * T**2 + Dcp[i] * T**3 + Ecp[i] * T**-2)
    return Cp

def calc_Cv():
    Cp = calc_Cp()
    Cv = []
    for i in range(len(comp)):
        Cv.append(Cp[i]/gamma[i])
    return Cv

def calc_Cv_mix():
    Cv = calc_Cv()
    Cv_mix = 0.0
    for i in range(len(comp)):
        Cv_mix = Cv_mix + comp[i]/100 * Cv[i] 
    return Cv_mix

def calc_Z():
    #this is an approximate calculation
    T_cr_mix_K = 0.0
    P_cr_mix_Pa = 0.0
    omega_mix = 0.0
    for i in range(len(comp)):
        T_cr_mix_K = T_cr_mix_K + comp[i]/100 * T_cr_K[i]     
        P_cr_mix_Pa = P_cr_mix_Pa + comp[i]/100 * P_cr_MPa[i]
        omega_mix = omega_mix + comp[i]/100 * omega[i]
    P_cr_mix_Pa = P_cr_mix_Pa * 1E6
    T_r = T/T_cr_mix_K
    P_r = P_a_Pa/P_cr_mix_Pa

    B0 = 0.083 - 0.422/(T_r**1.6)
    B1 = 0.139 - 0.172/(T_r**4.2)
    Z = 1 + (B0 + omega_mix * B1)*P_r/T_r
    T_r_lim = 0.686 + 0.439*P_r
    Check_T_r = "NO"
    if T_r > T_r_lim:
        Check_T_r = "OK"
    return Z, Check_T_r

def calc_dens():
    dens = MM_tot * P_a_Pa / (R * T * 1000)
    return dens

#data
P_env = 1.01 #in bar
T_C = -280.0 #inizialization for the control
P_barg = -1 #inizialization for the control
R = 8.3143
MM_H = 1.01
MM_O = 15.9994
MM_C = 12.001
MM_N = 14

order = ["H2","CO","CO2","CH4","N2","H2O","O2","C2H6","C3H8","C4H10","C5H12"]
comp = []
Acp = [6.62, 6.6, 10.34, 5.34, 6.5, 8.22, 8.27, 1.648, -0.966, 0.20372, -0.863]
Bcp = [0.00081, 0.0012, 0.00274, 0.0115, 0.001, 0.00015, 0.00026, 0.04142, 0.07279, 0.091307,0.0116]
Ccp = [0, 0, 0, 0, 0, 1.34E-6, 0, -1.5E-5, -3.8E-5, -4.6E-5, -5.92E-5]
Dcp = [0, 0, 0, 0, 0, 0, 0, 1.74E-9, 7.58E-9, 8.9E-9, 1.16E-8]
Ecp = [0, 0, -195500, 0, 0, 0, -187700, 0, 0, 0, 0]
gamma = [1.41, 1.4, 1.28, 1.32, 1.4, 1.33, 1.4, 1.18, 1.13, 1.18, 1.08]
T_cr_K = [41.25, 134.15, 304.25, 190.65, 126.05, 647.3, 154.35, 305.25, 369.95, 425, 470]
P_cr_MPa = [2.10756, 3.546, 7.397, 4.64068, 3.394388, 22.04832, 5.035853, 4.94466, 4.256, 3.8, 3.38]
omega = [-0.22, 0.049, 0.225, 0.008, 0.04, 0.344, 0.021, 0.098, 0.152, 0.193, 0.251]

#Auto calculation of MM
MM = []
MM.append(2*MM_H) #H2
MM.append(MM_C + MM_O) #CO
MM.append(MM_C + 2*MM_O) #CO2
MM.append(MM_C + 4*MM_H) #CH4
MM.append(2*MM_N) #N2
MM.append(MM_O + 2*MM_H) #H2O
MM.append(2*MM_O) #O2 
MM.append(2*MM_C + 6*MM_H) #C2H6
MM.append(3*MM_C + 8*MM_H) #C3H8
MM.append(4*MM_C + 10*MM_H) #C4H10
MM.append(5*MM_C + 12*MM_H) #C5H12

#Percentage
H2 = CO = CO2 = CH4 = N2 = H2O = O2 = C2H6 = C3H8 = C4H10 = C5H12 = 0.0
print("Insert the percentage of the components of the mix (if left empty, a zero will be assigned):")

for i in range(len(order)):
    print(str(order[i]) + " percentage: ")
    inp = input()
    if inp == "" or inp == " ":
        inp = 0.0
    else:
        inp = float(inp)
    comp.append(inp) 

print("Data inserted:")
print(f"Component\t%")
for i in range(len(order)):
    print(f"{order[i]}\t{comp[i]}")

#sum of percentage and check
sum_per = 0.0
for i in range(len(comp)):
    sum_per = sum_per + comp[i]

if round(sum_per,3) != 100.000 :
    print("Warning! The sum of the percentages is not 100.000")

while T_C < -273.15 :
    T_C = float(input("insert the Temperature (°C): "))
T = T_C + 273.15

while P_barg <= 0 :
    P_barg = float(input("insert the Pressure (barg): "))
P_bar = P_barg + P_env
P_a_Pa = P_bar * 1E5

MM_tot = MM_tot()
Cp_mix = calc_Cp_mix()
Cv_mix = calc_Cv_mix()
Z, Check_T_r = calc_Z()
dens = calc_dens()

#results
print("MM of the mixture is:", aprox(MM_tot))
print("Actual density of the mixture is:", aprox(dens), "kg/m^3")
print("Cp of the mixture is:", aprox(Cp_mix), "cal/molK")
print("Cv of the mixture is:", aprox(Cv_mix), "cal/molK")
print("Gamma of the mixture is:", Cp_mix/Cv_mix)
print("Approximate Compression of the mixture is:", aprox(Z))
print("And its control status is:", aprox(Check_T_r))
