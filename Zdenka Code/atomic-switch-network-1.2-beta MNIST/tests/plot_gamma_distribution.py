#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import scipy.special as sps
import numpy as np


number_of_wires = 2048
mean_len  = 30.0
std_len   = 10.0

gamma_shape = (mean_len/std_len)**2
gamma_scale = std_len**2 / mean_len

wire_lengths = np.random.gamma(gamma_shape, gamma_scale, number_of_wires)

count, bins, ignored = plt.hist(wire_lengths, 50, 
                              normed=True, 
                              histtype='stepfilled',
                              alpha = 0.75)


y = bins**(gamma_shape-1)*(np.exp(-bins/gamma_scale) / (sps.gamma(gamma_shape)*gamma_scale**gamma_shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.xlabel(r'Nanowire length [$\mu$m]')
plt.show()
