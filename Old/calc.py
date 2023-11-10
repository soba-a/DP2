import sympy as sp
import numpy as np
import math as m
from functions import *

# WORK IN N AND MM

# Known loads
Fp = sp.Matrix([-130.7, 0, 0]).T
Fb = sp.Matrix([42.8 * np.cos(np.radians(180 - 111.4)), 0, 42.8 * np.sin(np.radians(180-111.4))]).T
Fa = sp.Matrix([0, 10, 0]).T
Tp = sp.Matrix([0, 1548, 0]).T
Tb = sp.Matrix([0, -1548, 0]).T

# Loads we are trying to solve for
r1x, r1y, r1z, r2x, r2y, r2z = sp.symbols('r1x r1y r1z r2x r2y r2z')
R1 = sp.Matrix([r1x, r1y, r1z]).T
R2 = sp.Matrix([r2x, r2y, r2z]).T.subs(r2y, 0) # cause we assume Fa = +10 N

#### Solve for Reactions at bearings
# Fundamental problem we are solving: sumF = 0, sumM @ R1 = 0, sumM @ R2 = 0

# we know sum of forces = 0
sumF = R1 + R2 + Fp + Fb + Fa # we know this equals 0
# print("sumF: ", sumF)

# we know moments about R1 = 0
bearing_width = 0.5 * 25.4 # 1/2" to mm
# Calculating realistic case, at middle of bearing
MR2_at_r1 = sp.Matrix([0, 100 - bearing_width, 0]).T.cross(R2)
MFp_at_r1 = sp.Matrix([0, -50 - bearing_width/2, 0]).T.cross(Fp)
MFb_at_r1 = sp.Matrix([0, 115 - bearing_width/2, 0]).T.cross(Fb)
MFa_at_r1 = sp.Matrix([0, 115 - bearing_width/2, -203.2]).T.cross(Fa)

sumM_at_r1 = MR2_at_r1 + MFp_at_r1 + MFb_at_r1 + MFa_at_r1 + Tp + Tb # we know this equals 0

print("========================== COPY FOR \"SOLVE FOR REACTIONS\" ==========================")
print("MR2_at_r1: ", MR2_at_r1)
print("MFp_at_r1: ", MFp_at_r1)
print("MFb_at_r1: ", MFb_at_r1)
print("MFa_at_r1: ", MFa_at_r1)
print("sumM_at_r1: ", sumM_at_r1)
print("========================== END ==========================")

# create systems of equations!
syseq = sp.Matrix([sumF, sumM_at_r1])

# solve for unknowns
sol_loads = sp.solve(syseq, [r1x, r1y, r1z, r2x, r2z])


# Replace with values
R1 = R1.subs(sol_loads)
R2 = R2.subs(sol_loads)
MR2_at_r1 = MR2_at_r1.subs(sol_loads)
MFp_at_r1 = MFp_at_r1.subs(sol_loads)
MFb_at_r1 = MFb_at_r1.subs(sol_loads)
MFa_at_r1 = MFa_at_r1.subs(sol_loads)
print("========================== COPY FOR \"SOLVE THE SYSTEM\" ==========================")
print("sumF: ", sumF)
print("sumM_at_r1: ", sumM_at_r1)
print("Load Solutions:", sol_loads)
print("R1: ", R1)
print("R2: ", R2)
print("========================== END ==========================")

# print("=========================== Known Loads & Values ===========================")
# print("R1: ", R1)
# print("R2: ", R2)
# print("Fp: ", Fp)
# print("Fb: ", Fb)
# print("Fa: ", Fa)
# print("MR2_at_r1: ", MR2_at_r1)
# print("MFp_at_r1: ", MFp_at_r1)
# print("MFb_at_r1: ", MFb_at_r1)
# print("MFa_at_r1: ", MFa_at_r1)
# print("Tp: ", Tp)
# print("Tb: ", Tb)
# print("============================================================================")

## Welcome to the land of shaft design! We love designing shafts.

# Properties that do not change
S_ut = 751.529 # MPa
Se_prime = S_ut / 2 # MPa
S_y = 406.791 # MPa
n_shaft = 1.5

# Retaining Ring Groove is probably max stress (methinks)
y = sp.symbols('y')
rV = sp.Matrix([0, y, 0]).T
eq_loc_of_int = Tp + rV.cross(Fp)
Ma = (eq_loc_of_int.subs(y, 50)[0]**2 + eq_loc_of_int.subs(y, 50)[2]**2)**0.5 # assuming right at bearing edge
Tm = eq_loc_of_int.subs(y, 50)[1] # assuming right at bearing edge
print("Ma at Retaining Ring:", Ma)
print("Tm at Retaining Ring:", Tm)

### First Pass w/ Shigley Assumptions
kf_prime = 5 
kfs_prime = 3
k_ax = 5

# S_e & Marin Factors: see notes for reasoning
S_e_fp = calc_Se(Se_prime, ka=1, kb=0.9, kc=1, kd=1, ke=0.814, kf=1)

d_fp = ((16 * n_shaft / np.pi)*(2 * kf_prime * Ma * S_e_fp**(-1) + 3**(0.5) * kfs_prime * Tm * S_ut**(-1)))**(1/3)
print("First Pass Diameter: ", d_fp / 25.4, "inches")

### Second Pass, Calculate CORRECT values
d_2p = 0.5 * 25.4 # shaft dia, mm
print("Second Pass Diameter:", d_2p / 25.4, "inches")

# S_e & Marin Factors: see notes for reasoning. kb is different!
S_e_2 = calc_Se(Se_prime, ka=1, kb=calc_Kb(d_2p), kc=1, kd=1, ke=0.814, kf=1) # with correct kb value
print("Se:", S_e_2)

# Groove dim for 0.5" retaining ring :
groove_diameter = 0.468 * 25.4 # mm
groove_width = 0.039 * 25.4
groove_depth = 0.016 * 25.4

# Find Kt AND Kts IN GRAPHS! HERE ARE THE VALUES YOU NEED TO LOOK AT:
a = groove_width
t = groove_depth # = (d_2p - groove_diameter)
r = groove_width / 10

a_t = a/t
r_t = r/t
print("a/t: ", a_t, "r/t: ", r_t)
# go with r/t = 0.2, because that assumes groove radius is sharper than recommended

# Find Kt & Kts based on a/t & r/t. Shigley A-15-16 & 17
# assume Kt = Kf for now

# did you update your Kt and Kts?
Kt = 4.1 # Kt
Kts = 2.5 # Kfs 
# Kf = calc_K_f(Kt, d_2p, S_ut)
# Kfs = calc_K_fs(Kts, d_2p, S_ut)

# Fatigue Yield check! @ bearing edge
sig_a_p = sig_a_prime(Kt, Ma, d_2p)
sig_m_p = sig_m_prime(Kts, Tm, d_2p)
n = goodman_safety(sig_a=sig_a_p, sig_m=sig_m_p, S_e=S_e_2, S_ut=S_ut)
print("Goodman Safety Factor: ", n)
n = asme_safety(sig_a=sig_a_p, sig_m=sig_m_p, S_e=S_e_2, S_y=S_y)
print("ASME Safety Factor: ", n)
n = yield_check(sig_a=sig_a_p, sig_m=sig_m_p, S_y=S_y)
print("Yield Check: ", n)
n = yield_check_vonmise(M=Ma, Kf=Kt, T=Tm, d=d_2p, Kfs=Kts, S_y=S_y)
print("Yield Check Von Mise:", n)
# check static yield check
n = static_check_groove(F=Fa[1], M=Ma, Kt=Kt, T=Tm, d=d_2p, Kts=Kts, S_y=S_y)
print("Static Yield Check: ", n)

### Third Pass!

d_3p = 3/8 * 25.4
print("Third Pass Diameter:", d_3p / 25.4, "inches")

# S_e & Marin Factors: see notes for reasoning. kb is different!
S_e_3 = calc_Se(Se_prime, ka=1, kb=calc_Kb(d_3p), kc=1, kd=1, ke=0.814, kf=1) # with correct kb value
print("Se:", S_e_3)

# Groove dim for 3/8" retaining ring :
groove_diameter = 0.352 * 25.4 # mm
groove_width = 0.029 * 25.4
groove_depth = 0.0115 * 25.4

# Find Kt AND Kts IN GRAPHS! HERE ARE THE VALUES YOU NEED TO LOOK AT:
a = groove_width
t = groove_depth # = (d_2p - groove_diameter)
r = groove_width / 10

a_t = a/t
r_t = r/t
print("a/t: ", a_t, "r/t: ", r_t)
# go with r/t = 0.2, because that assumes groove radius is sharper than recommended

# Find Kt & Kts based on a/t & r/t. Shigley A-15-16 & 17
# assume Kt = Kf for now

# did you update your Kt and Kts?
Kt = 4.02 # Kt
Kts = 2.5 # Kfs 
# Kt = calc_K_f(Kt, d_2p, S_ut)
# Kts = calc_K_fs(Kts, d_2p, S_ut)

# Fatigue Yield check! @ bearing edge
sig_a_p = sig_a_prime(Kt, Ma, d_3p)
sig_m_p = sig_m_prime(Kts, Tm, d_3p)
n = goodman_safety(sig_a=sig_a_p, sig_m=sig_m_p, S_e=S_e_3, S_ut=S_ut)
print("Goodman Safety Factor: ", n)
n = asme_safety(sig_a=sig_a_p, sig_m=sig_m_p, S_e=S_e_3, S_y=S_y)
print("ASME Safety Factor: ", n)
n = yield_check(sig_a=sig_a_p, sig_m=sig_m_p, S_y=S_y)
print("Yield Check: ", n)
n = yield_check_vonmise(M=Ma, Kf=Kt, T=Tm, d=d_3p, Kfs=Kts, S_y=S_y)
print("Yield Check Von Mise:", n)
# check static yield check
n = static_check_groove(F=Fa[1], M=Ma, Kt=Kt, T=Tm, d=d_3p, Kts=Kts, S_y=S_y)
print("Static Yield Check: ", n)
von_mises = von_mises_groove(F=Fa[1], M=Ma, Kt=Kt, T=Tm, d=d_3p, Kts=Kts)
print("Von Mises at Groove:", von_mises, "MPa")

## COMPONENT DESIGN
# Key
S_ut_key = 515 #Mpa
S_y_key = S_ut_key * 0.8 # def

w = 1/8 * 25.4
h = w
# Shear Failure
r = d_2p / 2
F = Tm/r
n = 1.2
print("Shear Force on Key:", F)
L = 2 * F / w * n / S_y_key
print("minimum key Length:", L / 25.4)
L2 = 2 * F / h / S_y_key * n
print("mimimum key length2:", L2 / 25.4)
# key length = 1/4 "

# Set Screw
# no math needed for now...

# Pin
d_pin = 3/32 * 25.4 # mm
Sy_pin = (1400) / (2 * (np.pi * (3/32)**2 / 4)) # psi
Sy_pin = Sy_pin / 145.037738 # MPa
print(Sy_pin, "MPa")
n_pin = Sy_pin / 2 * 2 * np.pi * d_pin**2 / F / 4
print("Safety Factor Pin:", n_pin)

# Retaining Ring

