# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 2-6: Relations among thermal design conditions in a PWR


# Given power rates
q1_melt = 70    # kW/m
q1_avg = 17.86  # KW/m

# Peaking factors
RADIAL_FLUX = 1.55
AXIAL_FLUX = 1.70
ENGR_UNCERTAINTY = 1.05
OVERPOWER = 1.15
total_peaking = RADIAL_FLUX*AXIAL_FLUX*ENGR_UNCERTAINTY*OVERPOWER

q1_max = total_peaking*q1_avg

print("Total peaking:", round(total_peaking, 3))
print("q' (max):     6", round(q1_max, 3), "kW/m")
print("q' (melt) / q' (max):", round(q1_melt/q1_max, 2))


