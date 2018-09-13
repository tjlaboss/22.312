# 22.312, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Problem 3-7: Effect of continuous refueling on decay heat

from scipy.integrate import quad as integrate

Q = 3000
TAU_18 = 18*30.4*24*3600     # 18 months, in seconds: newer fuel
TAU_36 = 2*TAU_18            # 3 years, in seconds: oldest fuel

# Equation 3.70d
P1 = lambda t, tau: 0.066*(t**(-0.2) - (t + tau)**(-0.2))   # 't' is time after shutdown
# Modified for online refueling: (t + TAU_s) is always TAU_s
P2 = lambda t, tau: 0.066*(t**(-0.2) - (tau)**(-0.2))

# times
t_s = (60, 3600, 3600*24, 3600*24*30.4, 3600*24*365)
labels = ("1 min:  ", "1 hour: ", "1 day:  ", "1 month: ", "1 year:  ")

for i in range(len(t_s)):
	t = t_s[i]
	l = labels[i]
	# Case 1 (batch refueling): Average the 18 month old and the 3 year old fuel
	case1 = round(Q*(P1(t, TAU_18) + P1(t, TAU_36))/2, 2)
	
	# Case 2 (online refueling): Average over a differential dtau
	# Integral[ P(t, tau) dtau ] / Integral[ dtau ]   {over (0, TAU_36)}
	integrand = lambda tau: P1(t, tau)
	num = integrate(integrand, 0, TAU_36)[0]
	denom = TAU_36  # dtau over (0, TAU_36)
	case2 = round(Q*num/denom, 2)
	
	print(l, "Case 1", case1, "MWt; \tCase 2", case2, "MWt")
