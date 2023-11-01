import sympy as sp
import numpy as np
import math as m

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
# Calculating worst case, at inside bearing edge
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

# ## Welcome to the land of shaft design! We love designing shafts.
# # First Pass
# # Retaining Ring
# kf_prime = 5
# kfs_prime = 3
# k_ax = 5
# n_design = 1.5

# d_fp = ((16 * n_design / np.pi)(2 * kf_prime * Ma / S_e + 3**0.5 * kfs_prime * Tm / S_ut))**1/3


# shaft_d = 1/2 * 25.4 # in to mm 
# Se_prime = 337.5 # MPa
# # Marin factors, see notes for reasoning
# k_a = 1
# k_b = 0.879 * shaft_d **(-0.107)
# k_c = 1
# k_d = 1
# k_e = 0.814 # 99%
# k_f = 0.75 #???

# Se = Se_prime * k_a * k_b * k_c * k_d * k_e * k_f
# print("Se: ", Se)

