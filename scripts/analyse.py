#!/bin/python
import os
import re
import csv
import datetime
from dateutil.parser import parse
import xml.etree.cElementTree as et
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d
import tkinter as tk
from tkinter import filedialog

version_number = "0.1.0"


class SpeFile(object):
    """.spe file object.

    Attributes:
        xdim: 64 bit integer.
            Length of the x dimension of the image.
        ydim: 64 bit integer.
            Length of the y dimension of the image.
        wavelength: 1d numpy array.
            Wavelength dimension of the spectra.
        image: 2d numpy array.
            Image intensity values.
        footer: cElementTree object.
            XML file footer.
        time: datetime object.
            When the image file was created.
        exposure: 64 bit integer.
            Duration of exposure in milliseconds.
        roi_spectra: 1d numpy array.
            Mean intensity of the ROI for each wavelength.
        background_spectra: 1d numpy array.
            Mean intensity of the background for each wavelength.
        background_corrected_spectra: 1d numpy array.
            Mean intensity of the ROI for each wavelength with the
            background subtracted.
        intensities: Dictionary {Wavelength:Intensity}
            Variables stored as 64 bit floats.
        background_corrected_intensities: Dictionary
        {Wavelength:Background corrected intensity}
            Variables stored as 64 bit floats.
    """

    def __init__(self, filename):
        """Initialises SpeFile with its filename and proceeds to load
        the image and metadata.

        Arguments:
            filename: string
                Valid path to an .spe file.
        """
        self.filename = os.path.abspath(filename)
        self.xdim = None
        self.ydim = None
        self.image = None
        self.footer = None
        self.wavelength = None
        self.time = None
        self.exposure = None
        self.cumulative_exposure = None
        self.roi_spectra = None
        self.background_spectra = None
        self.background_corrected_spectra = None
        self.intensities = None
        self.bleaching_corrected_spectra = None
        self.bleaching_corrected_intensities = None
        self.bleaching_corrected_image = None
        self.concentration = None
        self.response = None

        self._load()

    def _load(self):
        """Loads image and relevant metadata from the file SpeObject was
        initialised with.
        """
        with open(self.filename, 'rb') as fid:
            self.xdim = np.int64(binary_to_array(fid, 42, 1, np.uint16)[0])
            self.ydim = np.int64(binary_to_array(fid, 656, 1, np.uint16)[0])
            self.image = np.float64(binary_to_array(fid, 4100, self.xdim * self.ydim, np.uint16)).reshape((self.ydim, self.xdim))
            self.footer = binary_to_xml(fid, binary_to_array(fid, 678, 1, np.uint64)[0])
            ns0 = '{http://www.princetoninstruments.com/spe/2009}'
            ns1 = '{http://www.princetoninstruments.com/experiment/2009}'
            self.time = parse(self.footer.find('.//{0}FileInformation'.format(ns0)).get('created'))
            self.exposure = np.int64(self.footer.find('.//{0}ExposureTime'.format(ns1)).text)
            self.wavelength = np.array([np.float64(s) for s in self.footer.find('.//{0}Wavelength'.format(ns0)).text.split(',')])

    def correct_background(self, roi, background):
        """Corrects ROI for background fluorescence.

        Arguments:
            roi: tuple or list of two integers.
                Provides the y range of the chose ROI.
            background: tuple or list of two integers.
                Provides the y range of the chose background region.
        """
        self.roi_spectra = np.mean(self.image[min(roi):max(roi), :], axis=0)
        self.background_spectra = np.mean(self.image[min(background):max(background), :], axis=0)
        # background_mask = np.ones((self.ydim, self.xdim), dtype=bool)
        # background_mask[min(roi):max(roi), :] = False
        # self.background_spectra = np.median(self.image[background_mask], axis=0)
        self.background_corrected_spectra = self.roi_spectra - self.background_spectra
        self.corrected_image = self.image[min(roi):max(roi), :] - self.background_spectra
        self.summed_spectra = np.sum(self.image[min(roi):max(roi), :] - self.background_spectra, axis=0)

    def calculate_intensities(self, peaks):
        self.intensities = dict()
        self.summed_intensities = dict()
        if not isinstance(peaks, list):
            peaks = [peaks]
        for index in peaks:
            wavelength = np.mean(self.wavelength[index - 15:index + 15])
            intensity = np.mean(self.background_corrected_spectra[index - 15:index + 15])
            self.intensities[wavelength] = intensity
            intensity = np.mean(self.summed_spectra[index - 15:index + 15])
            self.summed_intensities[wavelength] = intensity

    def correct_bleaching(self, wavelength, factor):
        if self.bleaching_corrected_intensities is None:
            self.bleaching_corrected_intensities = dict()
        if self.bleaching_corrected_spectra is None:
            self.bleaching_corrected_spectra = dict()
        self.bleaching_corrected_spectra[wavelength] = self.background_corrected_spectra / factor
        self.bleaching_corrected_intensities[wavelength] = self.intensities[wavelength] / factor
        if wavelength < 480:
            self.bleaching_corrected_image = self.corrected_image / factor


class RoiPicker(object):
    """Pick ROIs in a 2d numpy array.

    Given a 2d numpy array, this object provides an interactive pyplot
    based GUI to select either horizontal or vertical ROIs.

    Attributes:
        indices: list of two integers.
            The indices of the two lines chosen; i.e. a slice of the
            array between these points gives the chosen ROI.
        axis: 8 bit integer.
            The axis of the ROI, provided on initialisation. 0
            represents a horizontal ROI, 1 a vertical.
    """

    def __init__(self, image=None, x=None, y=None, axis=None):
        """Initialise a RoiPicker object with an image and a chosen axis.

        Arguments:
            image: 2d numpy array.
                Array to choose a ROI from.
            axis: {0,1}.
                Axis along which to select a ROI. 0 represents a
                horizontal ROI, 1 a vertical.
        """
        self.indices = []
        self.axis = np.int8(axis)
        self._figure = plt.figure(figsize=(14, 4))
        self._axes = self._figure.add_subplot(111)
        if image is not None:
            self._axes.imshow(image)
        elif x is not None and y is not None:
            self._axes.plot(x, y)
        if self.axis == 0:
            self._x = self._axes.get_xlim()
        elif self.axis == 1:
            self._y = self._axes.get_ylim()
        self._connection = self._figure.canvas.mpl_connect('button_press_event', self._button_press_callback)
        plt.show(block=True)

    def _button_press_callback(self, event):
        if event.inaxes:
            if self.axis == 0:
                self._y = (event.ydata, event.ydata)
                self.indices.append(np.int64(event.ydata))
            elif self.axis == 1:
                self._x = (event.xdata, event.xdata)
                self.indices.append(np.int64(event.xdata))
            self._axes.plot(self._x, self._y, color="red", linewidth=1.0)
            self._figure.canvas.draw()
            if len(self.indices) == 2:
                self._figure.canvas.mpl_disconnect(self._connection)
                plt.close()


def load_files(path):
    files = []
    if isinstance(path, str) is True:
        if os.path.isfile(path) is True and re.search('.spe', path):
            files.append(SpeFile(path))
        elif os.path.isdir(path) is True:
            for subpath in os.listdir(path):
                if re.search('.spe', subpath):
                    files.append(SpeFile(os.path.join(path, subpath)))
    elif isinstance(path, tuple) is True or isinstance(path, list) is True:
        for entry in path:
            if re.search('.spe', entry):
                files.append(SpeFile(entry))
    files.sort(key=lambda f: f.time)
    exposure = 0
    for file in files:
        file.cumulative_exposure = exposure
        exposure += file.exposure
    return files


def subtract_background(files, index):
    roi = RoiPicker(image=files[index].image, axis=0)
    background = RoiPicker(image=files[index].image, axis=0)
    for file in files:
        file.correct_background(roi.indices, background.indices)
    return files


def find_peaks(file, range=False, smooth=True, mph=20, mpd=30):
    """Detect peaks in the background corrected spectra.
    Arguments:
        smooth: boolean, optional (default = True)
            If true, smooth the spectra with a savgol filter to minimise
            false positives. Recommended.
        mph: {None, int}, optional (default = 20)
            Only returns peaks greater than this height.
        mpd: {None, int}, optional (default = 30)
            Sets a minimum distance between peaks to return.
    """
    if range is True:
        roi = RoiPicker(x=file.wavelength, y=file.background_corrected_spectra, axis=1)
        spectra = file.background_corrected_spectra[min(roi.indices):max(roi.indices)]
    else:
        spectra = file.background_corrected_spectra
    if smooth is True:
        array = savgol_filter(spectra, 31, 3)
    else:
        array = spectra
    dx = array[1:] - array[:-1]
    peaks = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    if mph is not None:
        peaks = peaks[array[peaks] >= mph]
    if mpd is not None and len(peaks) > 1:
        peaks = peaks[np.argsort(array[peaks])][::-1]  # sorts peaks by magnitude
        idel = np.zeros(peaks.size, dtype=bool)
        for i in range(peaks.size):  # checks for peaks closer than the mpd
            if not idel[i]:
                idel = idel | (peaks >= peaks[i] - mpd) & (peaks <= peaks[i] + mpd)
                idel[i] = 0
        peaks = np.sort(peaks[~idel])  # resorts peaks and removes those that are too close
    return peaks


def bleaching(files, function):
    values = dict()
    exposure = np.array([file.cumulative_exposure / 1000 for file in files])
    for wavelength in file.intensities.keys():
        intensities = np.array([file.intensities[wavelength] for file in files])
        normalised_intensities = np.array([intensity / intensities[0] for intensity in intensities])
        popt, pcov = curve_fit(function, exposure, normalised_intensities, p0=(1, 0.05), maxfev=2000)
        plt.plot(exposure, normalised_intensities, 'bo', exposure, function(exposure, *popt), 'r-')
        plt.title(str(wavelength))
        plt.show(block=True)
        values[wavelength] = popt
    return values


def show_images(files, *args, **kwargs):
    for file in files:
        fig = plt.figure(figsize=(14, 4))
        axes = fig.add_subplot(111)
        axes.imshow(file.corrected_image, *args, **kwargs)
        axes.set_title(file.filename)
        plt.show()


def binary_to_array(file, offset, size, dtype):
    file.seek(offset)
    return np.fromfile(file, dtype, size)


def binary_to_xml(file, offset):
    file.seek(offset)
    return et.parse(file)


def get_index_from_wavelength(file, wavelength):
    return min(range(len(file.wavelength)), key=lambda i: abs(file.wavelength[i] - wavelength))


def bleaching_func(x, a, b):
    return a * np.exp(-b * x) + (1 - a)


if __name__ == "__main__":
    plt.ion()
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilenames()
    root.destroy()
    testexp = load_files(path)
    testexp_bck = subtract_background(testexp, 0)
    for file in testexp_bck:
        plt.plot(file.wavelength, file.background_corrected_spectra)
    plt.show(block=True)
    fig, axes = plt.subplots(len(testexp_bck), 1, sharex=True, sharey=True)
    for i, file in enumerate(testexp_bck, 0):
        image = file.corrected_image[:, 400:1000]
        if i == 0:
            lowest_intensity = min(image.flatten())
            highest_intensity = max(image.flatten())
        wavelength = file.wavelength[400:1000]
        axes[i].imshow(
            image,
            extent=[min(wavelength), max(wavelength), 0, 20],
            cmap="inferno",
            vmin=lowest_intensity,
            vmax=highest_intensity
        )
    plt.show(block=True)
    peaks = [get_index_from_wavelength(testexp_bck[0], i) for i in [472, 508, 561]]
    for index, file in enumerate(testexp_bck):
        file.calculate_intensities(peaks)
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilenames()
    root.destroy()
    bl_files = []
    for i in testexp_bck:
        if i.filename in path:
            bl_files.append(i)
    xdict = bleaching(bl_files, bleaching_func)
    intensities = []
    corrected_intensities = []
    for file in testexp_bck:
        for wavelength in xdict.keys():
            popt = xdict[wavelength]
            norm = bleaching_func(file.cumulative_exposure / 1000, *popt)
            file.correct_bleaching(wavelength, norm)
    fig, axes = plt.subplots(len(testexp_bck), 1, sharex=True, sharey=True)
    for i, file in enumerate(testexp_bck, 0):
        image = file.bleaching_corrected_image[:, 400:1100]
        if i == 0:
            lowest_intensity = min(image.flatten())
            highest_intensity = max(image.flatten())
        wavelength = file.wavelength[400:1100]
        axes[i].imshow(
            image,
            extent=[min(wavelength), max(wavelength), 0, 20],
            cmap="inferno",
            vmin=lowest_intensity,
            vmax=highest_intensity
        )
    plt.show(block=True)
    fmax = dict()
    for wavelength in testexp_bck[0].bleaching_corrected_intensities.keys():
        fmax[wavelength] = testexp_bck[0].bleaching_corrected_intensities[wavelength]
    for file in testexp_bck:
        conc = input("Enter concentration (M) for file {0}:".format(file.filename))
        if conc == "":
            file.concentration = None
        else:
            file.concentration = float(conc)
        file.response = dict()
        for wavelength in file.bleaching_corrected_intensities.keys():
            file.response[wavelength] = (1 - (file.bleaching_corrected_intensities[wavelength] / fmax[wavelength]))
    conc = [file.concentration for file in testexp_bck]
    responses = []
    for wavelength in testexp_bck[0].bleaching_corrected_intensities.keys():
        for file in testexp_bck:
            if wavelength < 480:
                responses.append(file.response[wavelength])
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for wavelength in testexp_bck[0].bleaching_corrected_spectra.keys():
        for file in testexp_bck:
            if file.concentration > 0 and wavelength < 480:
                ax1.plot(file.wavelength, file.bleaching_corrected_spectra[wavelength], label=str(file.concentration))
    ax1.legend()
    ax2.plot(responses, 'bo')
    plt.show(block=True)
    construct = input("Enter name of construct:")
    date = testexp_bck[0].time.date().isoformat()
    root = tk.Tk()
    root.withdraw()
    path = filedialog.asksaveasfilename(initialfile=construct + "_" + date)
    path_intensities = path + "_intensities.csv"
    path_spectra = path + "_spectra.csv"
    path_images = path + "_images.csv"
    root.destroy()

    with open(path_intensities, 'w', newline='') as f:
        fieldnames = ['construct', 'date', 'cumulative_exposure', 'concentration', 'dye', 'response', 'raw_intensity', 'corrected_intensity']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for file in testexp_bck:
            for dye in file.bleaching_corrected_intensities.keys():
                results = dict()
                results["dye"] = dye
                results["construct"] = construct
                results["date"] = date
                results["cumulative_exposure"] = file.cumulative_exposure
                results["concentration"] = file.concentration
                if file.response is not None:
                    results["response"] = file.response[dye]
                else:
                    results["response"] = None
                results["raw_intensity"] = file.intensities[dye]
                results["corrected_intensity"] = file.bleaching_corrected_intensities[dye]
                writer.writerow(results)

    with open(path_spectra, 'w', newline='') as f:
        fieldnames = ['construct', 'date', 'cumulative_exposure', 'concentration', 'wavelength', 'dye', 'raw_spectra', 'background_corrected_spectra', 'bleaching_corrected_spectra']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for file in testexp_bck:
            for dye in file.bleaching_corrected_intensities.keys():
                for index, wavelength in enumerate(file.wavelength):
                    results = dict()
                    results["dye"] = dye
                    results["construct"] = construct
                    results["date"] = date
                    results["cumulative_exposure"] = file.cumulative_exposure
                    results["concentration"] = file.concentration
                    results["wavelength"] = wavelength
                    results["raw_spectra"] = file.roi_spectra[index]
                    results["background_corrected_spectra"] = file.background_corrected_spectra[index]
                    results["bleaching_corrected_spectra"] = file.bleaching_corrected_spectra[dye][index]
                    writer.writerow(results)

    with open(path_images, 'w', newline='') as f:
        fieldnames = ['construct', 'date', 'cumulative_exposure', 'concentration', 'wavelength', 'y', 'bleaching_corrected_intensity']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for file in testexp_bck:
            for x_index, wavelength in enumerate(file.wavelength):
                for y_index, intensity in enumerate(file.bleaching_corrected_image[:, x_index]):
                    results = dict()
                    results["construct"] = construct
                    results["date"] = date
                    results["cumulative_exposure"] = file.cumulative_exposure
                    results["concentration"] = file.concentration
                    results["wavelength"] = wavelength
                    results["y"] = y_index
                    results["bleaching_corrected_intensity"] = intensity
                    writer.writerow(results)
