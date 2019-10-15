#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.exposure import rescale_intensity
from skimage import img_as_float
from skimage.transform import rotate

wtconstructs = img_as_float(rescale_intensity(io.imread("data/180612_wt_western.tif", plugin="tifffile", as_gray=True), out_range="float64"))
gfpconstructs = img_as_float(rescale_intensity(io.imread("data/180612_gfp_western.tif", plugin="tifffile", as_gray=True), out_range="float64"))

wtconstructs = rotate(wtconstructs, 357.5)[90:340, 100:300]
gfpconstructs = rotate(gfpconstructs, 357.5)[100:350, 110:310]

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(wtconstructs, cmap="Greys")
ax2.imshow(gfpconstructs, cmap="Greys")
plt.show(block=True)


def plot_profile(array, i1, i2):
    block = np.array(array[:, i1:i2])
    profile = np.mean(block, axis=1)
    return(profile)


background = np.mean(wtconstructs[:, 170:])
wt = rescale_intensity(wtconstructs[:, :40] - background, out_range=(0, 1))
w311anapsur = rescale_intensity(wtconstructs[:, 42:80] - background, out_range=(0, 1))
w311anap = rescale_intensity(wtconstructs[:, 80:120] - background, out_range=(0, 1))
w311 = rescale_intensity(wtconstructs[:, 120:160] - background, out_range=(0, 1))

backgroundgfp = np.mean(gfpconstructs[:, 170:])
wtgfp = rescale_intensity(gfpconstructs[:, 4:47] - background, out_range=(0, 1))
w311anapsurgfp = rescale_intensity(gfpconstructs[:, 48:84] - background, out_range=(0, 1))
w311anapgfp = rescale_intensity(gfpconstructs[:, 84:124] - background, out_range=(0, 1))
w311gfp = rescale_intensity(gfpconstructs[:, 120:160] - background, out_range=(0, 1))

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
ax1.imshow(wt, cmap="Greys")
ax2.imshow(w311anapsur, cmap="Greys")
ax3.imshow(w311anap, cmap="Greys")
ax4.imshow(w311, cmap="Greys")
plt.show(block=True)

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(1, 8, sharey="row")
ax1.imshow(wtgfp, cmap="Greys", aspect="auto")
ax3.imshow(w311anapsurgfp, cmap="Greys", aspect="auto")
ax5.imshow(w311anapgfp, cmap="Greys", aspect="auto")
ax7.imshow(w311gfp, cmap="Greys", aspect="auto")

y1 = np.mean(wtgfp, axis=1)
y2 = np.mean(w311anapsurgfp, axis=1)
y3 = np.mean(w311anapgfp, axis=1)
y4 = np.mean(w311gfp, axis=1)
x = np.arange(0, len(y1), 1)
a = x > 50
b = x < 150
c = x > 150
d = x < 250
full_length_mask = np.all([a, b], axis=0)
truncated_mask = np.all([c, d], axis=0)

full_length = np.trapz(y1[50:150])
truncated = np.trapz(y1[150:250])
full_length_frac = (full_length / (full_length + truncated)) * 100
truncated_frac = (truncated / (full_length + truncated)) * 100

ax2.plot(y1, x, 'k-')
ax2.invert_yaxis()
ax2.fill_between(y1, 0, x, where=full_length_mask, facecolor='blue')
ax2.annotate("{0}%".format(round(full_length_frac, 1)), xy=(0.3, 50))
ax2.fill_between(y1, 0, x, where=truncated_mask, facecolor='orange')
ax2.annotate("{0}%".format(round(truncated_frac, 1)), xy=(0.3, 200))

full_length = np.trapz(y2[50:150])
truncated = np.trapz(y2[150:250])
full_length_frac = (full_length / (full_length + truncated)) * 100
truncated_frac = (truncated / (full_length + truncated)) * 100

ax4.plot(y2, x, 'k-')
ax4.invert_yaxis()
ax4.fill_between(y2, 0, x, where=full_length_mask, facecolor='blue')
ax4.annotate("{0}%".format(round(full_length_frac, 1)), xy=(0.3, 50))
ax4.fill_between(y2, 0, x, where=truncated_mask, facecolor='orange')
ax4.annotate("{0}%".format(round(truncated_frac, 1)), xy=(0.3, 200))

full_length = np.trapz(y3[50:150])
truncated = np.trapz(y3[150:250])
full_length_frac = (full_length / (full_length + truncated)) * 100
truncated_frac = (truncated / (full_length + truncated)) * 100

ax6.plot(y3, x, 'k-')
ax6.invert_yaxis()
ax6.fill_between(y3, 0, x, where=full_length_mask, interpolate=True, facecolor='blue')
ax6.annotate("{0}%".format(round(full_length_frac, 1)), xy=(0.3, 50))
ax6.fill_between(y3, 0, x, where=truncated_mask, interpolate=True, facecolor='orange')
ax6.annotate("{0}%".format(round(truncated_frac, 1)), xy=(0.3, 200))

full_length = np.trapz(y4[50:150])
truncated = np.trapz(y4[150:250])
full_length_frac = (full_length / (full_length + truncated)) * 100
truncated_frac = (truncated / (full_length + truncated)) * 100

ax8.plot(y4, x, 'k-')
ax8.invert_yaxis()
ax8.fill_between(y4, 0, x, where=full_length_mask, interpolate=True, facecolor='blue')
ax8.annotate("{0}%".format(round(full_length_frac, 1)), xy=(0.3, 50))
ax8.fill_between(y4, 0, x, where=truncated_mask, interpolate=True, facecolor='orange')
ax8.annotate("{0}%".format(round(truncated_frac, 1)), xy=(0.3, 200))

plt.show(block=True)
