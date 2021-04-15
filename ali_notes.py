"""
Notes Re: Berend's Code
"""

# x is energy
# y is
# x and y are always what you would see on the detector (computer screen) as you measure things in the lab
# z is what you're measuring / scanning through

# ARPES detector is 2D, 1 energy axis, 1 angle axis (theta in Berend's code)
# then you always scan in
# z is often the other angle axis (SES phi -- related to voltages set in analyser), but could be temperature dependence,
# photon energy dependence...
# SES phi is simulating the angle of the sample (off perpendicular) to the slit. So at a given SES phi, you scan a
# slice of energy vs. phi along a given kx-ky. But changing SES phi lets us simulate changing the angle of the sample
# which puts us at a different location of kx-ky


# plt.plot can plot array as a function of indices, or as a function of the x-axis scaling
