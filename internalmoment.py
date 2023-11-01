from calc import *
import matplotlib.pyplot as plt

# # time to just make an internal moment graph
# y = sp.symbols('y')
# rV = sp.Matrix([0, y, 0]).T
# eq_0_50 = Tp + rV.cross(Fp)
# eq_0_50 = eq_0_50.norm()

# eq_50_150 = Tp + rV.cross(Fp) + (rV - sp.Matrix([0, 50, 0]).T).cross(R1)
# eq_50_150 = eq_50_150.norm()

# eq_150_165 = Tp + rV.cross(Fp) + (rV - sp.Matrix([0, 50, 0]).T).cross(R1) + (rV - sp.Matrix([0, 150, 0]).T).cross(R2)
# eq_150_165 = eq_150_165.norm()

# # Plotting the equations
# y_vals = np.linspace(0, 165, 1000)
# eq_0_50_vals = [eq_0_50.subs(y, y_val) for y_val in y_vals if y_val <= 50]
# eq_50_150_vals = [eq_50_150.subs(y, y_val) for y_val in y_vals if 50 < y_val <= 150]
# eq_150_165_vals = [eq_150_165.subs(y, y_val) for y_val in y_vals if y_val > 150]

# plt.plot(y_vals[:len(eq_0_50_vals)], eq_0_50_vals, label='0 <= y <= 50')
# plt.plot(y_vals[len(eq_0_50_vals):len(eq_0_50_vals)+len(eq_50_150_vals)], eq_50_150_vals, label='50 < y <= 150')
# plt.plot(y_vals[len(eq_0_50_vals)+len(eq_50_150_vals):], eq_150_165_vals, label='150 < y <= 165')

# plt.xlabel('y (mm)')
# plt.ylabel('Internal Moment (Nmm)')
# plt.title('Internal Moment vs. y (LEFT TO RIGHT)')
# plt.legend()
# plt.show()

# # now do the other side (verified correct)
# eq_0_15 = Tb + rV.cross(Fb) + (rV - sp.Matrix([0, 0, 203.2]).T).cross(Fa)
# eq_0_15 = eq_0_15.norm()

# eq_15_115 = Tb + rV.cross(Fb) + (rV - sp.Matrix([0, 0, 203.2]).T).cross(Fa) + (rV - sp.Matrix([0, 15, 0]).T).cross(R2)
# eq_15_115 = eq_15_115.norm()

# eq_115_165 = Tb + rV.cross(Fb) + (rV - sp.Matrix([0, 0, 203.2]).T).cross(Fa) + (rV - sp.Matrix([0, 15, 0]).T).cross(R2) + (rV - sp.Matrix([0, 115, 0]).T).cross(R1)
# eq_115_165 = eq_115_165.norm()

# y_vals = np.linspace(0, 165, 1000)
# eq_0_15_vals = [eq_0_15.subs(y, y_val) for y_val in y_vals if y_val <= 15]
# eq_15_115_vals = [eq_15_115.subs(y, y_val) for y_val in y_vals if 15 < y_val <= 115]
# eq_115_165_vals = [eq_115_165.subs(y, y_val) for y_val in y_vals if y_val > 115]

# plt.plot(y_vals[:len(eq_0_15_vals)], eq_0_15_vals, label='0 <= y <= 15')
# plt.plot(y_vals[len(eq_0_15_vals):len(eq_0_15_vals)+len(eq_15_115_vals)], eq_15_115_vals, label='15 < y <= 115')
# plt.plot(y_vals[len(eq_0_15_vals)+len(eq_15_115_vals):], eq_115_165_vals, label='115 < y <= 165')

# plt.xlabel('y (mm)')
# plt.ylabel('Internal Moment (Nmm)')
# plt.title('Internal Moment vs. y (RIGHT TO LEFT)')
# plt.legend()
# plt.show()


## you only have tota moments! split into bending moment and torque
# time to just make an internal moment graph
y = sp.symbols('y')
rV = sp.Matrix([0, y, 0]).T
eq_0_b1 = Tp + rV.cross(Fp)
eq_0_b1_bend = sp.sqrt(eq_0_b1[0]**2 + eq_0_b1[2]**2)
eq_0_b1_tors = eq_0_b1[1]

eq_b1_b2 = Tp + rV.cross(Fp) + (rV - sp.Matrix([0, 50 + bearing_width/2, 0]).T).cross(R1)
eq_b1_b2_bend = sp.sqrt(eq_b1_b2[0]**2 + eq_b1_b2[2]**2)
eq_b1_b2_tors = eq_b1_b2[1]

eq_b2_165 = Tp + rV.cross(Fp) + (rV - sp.Matrix([0, 50 + bearing_width/2, 0]).T).cross(R1) + (rV - sp.Matrix([0, 150 - bearing_width/2, 0]).T).cross(R2)
eq_b2_165_bend = sp.sqrt(eq_b2_165[0]**2 + eq_b2_165[2]**2)
eq_b2_165_tors = eq_b2_165[1]

# Plotting the equations
y_vals = np.linspace(0, 165, 1000)
eq_0_b1_bend_vals = [eq_0_b1_bend.subs(y, y_val) for y_val in y_vals if y_val <= 50+bearing_width/2]
eq_b1_b2_bend_vals = [eq_b1_b2_bend.subs(y, y_val) for y_val in y_vals if 50+bearing_width/2 < y_val <= 150-bearing_width/2]
eq_b2_165_bend_vals = [eq_b2_165_bend.subs(y, y_val) for y_val in y_vals if y_val > 150-bearing_width/2]
eq_0_b1_tors_vals = [eq_0_b1_tors.subs(y, y_val) for y_val in y_vals if y_val <= 50+bearing_width/2]
eq_b1_b2_tors_vals = [eq_b1_b2_tors.subs(y, y_val) for y_val in y_vals if 50+bearing_width/2 < y_val <= 150-bearing_width/2]
eq_b2_165_tors_vals = [eq_b2_165_tors.subs(y, y_val) for y_val in y_vals if y_val > 150-bearing_width/2]

fig, axs = plt.subplots(1, 2, figsize=(12, 12))

axs[0].plot(y_vals[:len(eq_0_b1_bend_vals)], eq_0_b1_bend_vals, label='0 <= y <= Bearing 1')
axs[0].plot(y_vals[len(eq_0_b1_bend_vals):len(eq_0_b1_bend_vals)+len(eq_b1_b2_bend_vals)], eq_b1_b2_bend_vals, label='Bearing 1 < y <= Bearing 2')
axs[0].plot(y_vals[len(eq_0_b1_bend_vals)+len(eq_b1_b2_bend_vals):], eq_b2_165_bend_vals, label='Bearing 2 < y <= 165')
axs[0].set_xlabel('y (mm)')
axs[0].set_ylabel('Bending Moment (Nmm)')
axs[0].set_title('Bending Moment vs. y (LEFT TO RIGHT)')
axs[0].legend()

axs[1].plot(y_vals[:len(eq_0_b1_tors_vals)], eq_0_b1_tors_vals, label='0 <= y <= Bearing 1')
axs[1].plot(y_vals[len(eq_0_b1_tors_vals):len(eq_0_b1_tors_vals)+len(eq_b1_b2_tors_vals)], eq_b1_b2_tors_vals, label='Bearing 1 < y <= Bearing 2')
axs[1].plot(y_vals[len(eq_0_b1_tors_vals)+len(eq_b1_b2_tors_vals):], eq_b2_165_tors_vals, label='Bearing 2 < y <= 165')
axs[1].set_xlabel('y (mm)')
axs[1].set_ylabel('Torque (Nmm)')
axs[1].set_title('Torque vs. y (LEFT TO RIGHT)')
axs[1].legend()

plt.show()