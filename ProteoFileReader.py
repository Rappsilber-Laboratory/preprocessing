# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 14:08:39 2014

@author: sven
"""


#==============================================================================
# Datastructures
#==============================================================================
import re
import numpy as np
import sys
import pyopenms as oms
import io
import os


def read_mgf(mgf_file):
    mgf_reader = io.open(mgf_file, "r")

    spectra = []
    peaks = []
    RT, pep_mz, pep_int, charge, scanID, ms2id = -1, -1, -1, -1, -1, -1
    title, detector, fragmethod = "", "", ""
    peak_re = re.compile(r'^([0-9.e+\-]+)\s([0-9.e+\-]+)')

    for line in mgf_reader:
        if len(line.strip()) == 0:
            continue
        if re.match(peak_re, line):
            mz_int_list = re.match(peak_re, line).groups()
            peaks.append([float(x) for x in mz_int_list])
        elif line.startswith("TITLE"):
            title = line.replace('TITLE=', '').strip()
            # scan_match = re.search("^TITLE=[^.]*.([0-9]+).", line)
            # ms2id_scan_match = re.search("ms2_scanId=([0-9]+)", line)
            # scanID = int(scan_match.groups()[0])
            # if ms2id_scan_match:
            #     ms2id = int(ms2id_scan_match.groups()[0])
            # else:
            #     ms2id = -1

        elif line.startswith("RTINSECONDS"):
            RT = float(re.search("RTINSECONDS=(.*)", line).groups()[0])

        elif line.startswith("PEPMASS"):
            precursor = re.search("PEPMASS=(.*)", line).groups()[0].split()
            pep_mz = float(precursor[0])
            try:
                pep_int = float(precursor[1])
            except:
                pep_int = -1.0

        elif line.startswith("CHARGE"):
            charge = float(re.search("CHARGE=(-?\d)", line).groups()[0])

        elif line.startswith("DETECTOR"):
            detector = re.search("DETECTOR=(.*)", line).groups()[0].strip()

        elif line.startswith("FRAGMETHOD"):
            fragmethod = re.search("FRAGMETHOD=(.*)", line).groups()[0].strip()

        elif "END IONS" in line:
            spectra.append(
                MS2_spectrum(title, RT, pep_mz, pep_int, charge, peaks, detector=detector, fragmethod=fragmethod)
            )
            peaks = []
            title, detector, fragmethod = "", "", ""
            RT, pep_mz, pep_int, charge, scanID, ms2id = -1, -1, -1, -1, -1, -1

    return spectra


def write_mgf(spectra, outfile):
    out_writer = open(os.path.join(outfile), "w")
    out_writer.write('MASS=Monoisotopic\n')
    for spectrum in spectra:
        title = spectrum.getTitle()
        # scan = re.search('scan=[0-9]*', title).group(0)[5:]
        # try:
        #     title = re.match('([A-Z])[0-9]{6}_[0-9]{2}.+?( )', title).group(0)[:-1]
        # except AttributeError:
        #     title = re.match('[0-9]{8}_[0-9]{2}.+?( )', title).group(0)[:-1]
        # title = '.'.join([title, scan, scan, str(int(spectrum.charge))])
        if 'ms2_scanId' in spectrum.getTitle():
            try:
                ms2_parent = re.search('ms2_scanId=.*scan=([0-9]+)', spectrum.getTitle()).groups()[0]
            except AttributeError:
                ms2_parent = 0
            title += ' ms2_scanId=%s' % ms2_parent
        stavrox_mgf = """
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
DETECTOR={}
FRAGMETHOD={}
{}
END IONS""".format(
                title,
                spectrum.getPrecursorMZ(),
                spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                int(spectrum.charge),
                spectrum.getRT(),
                spectrum.getDetector(),
                spectrum.getFragMethod(),
                "\n".join([f"{mz} {i}" for mz, i in spectrum.peaks if i > 0])
            )
        out_writer.write(stavrox_mgf)


def split_mgf_methods(mgf_in_file):
    ms2_spectra = read_mgf(mgf_in_file)

    methods = [
        "CID",
        "HCD",
        "ETD",
        "ETciD",
        "EThcD"
    ]

    for method in methods:
        split_spectra = [spectrum for spectrum in ms2_spectra if spectrum.getFragMethod() == method]

        if len(split_spectra) > 0:
            out_file_name = '%s_%s' % (method, os.path.split(mgf_in_file)[1])
            out_file_path = os.path.join(os.path.split(mgf_in_file)[0], out_file_name)

            write_mgf(split_spectra, out_file_path)


def mzMLReader(in_file):
    """
    One line wrapper for OpenMS mzML reading. Returns the "exp" of a file.
    
    Parameters:
    -----------------------
    in_file: str, 
              location of the mzML file.
    """
    file = oms.MzMLFile()
    exp = oms.MSExperiment()
    file.load(in_file, exp)
    return exp


class MS2_spectrum():
    """
    Class container for MS2 spectra.
    We need the following input parameters:
    title, RT, pepmass, pepint, charge, peaks, peakcharge=[]

    Parameters:
    -----------------------------------------
    title: str,
            title of the spectrum
    RT: float,
        retention time of the precursor
    pepmass: float,
              mass of the precursor
    charge: int,
             charge of the precursor
    peaks: [(float, float)],
           mass intensity
    peakcharge: arr,
                charge array for the peaks

    """
    def __init__(self, title, RT, pepmz, pepint, charge, peaks, peakcharge=[], fragmethod='', detector=''):
        self.title = title
        self.RT = RT
        self.pepmz = pepmz
        self.pepint = pepint
        self.charge = charge
        self.peaks = peaks
        self.peakcharge = peakcharge
        self.fragMethod = fragmethod
        self.detector = detector

    def getPrecursorMZ(self):
        """
        Returns the precursor mass
        """
        return self.pepmz

    def getPrecursorIntensity(self):
        """
        Returns the precursor intensity
        """
        return self.pepint

    def getRT(self):
        """
        Returns the precursor RT
        """
        return self.RT

    def getTitle(self):
        """
        Returns the precursor mass
        """
        return self.title

    def getPeaks(self):
        """
        Returns the spectrum peaks
        """
        return self.peaks

    def getMZ(self):
        """
        Returns the mz of the MS2
        """
        return self.peaks[:, 0]

    def getIntensities(self):
        """
        Returns the MS2 peak intensities
        """
        return self.peaks[:, 1]

    def getUnchargedMass(self):
        """
        Computs the uncharged mass of a fragment:
        uncharged_mass = (mz - hydrogen mass) * z
        """
        return (self.pepmz - 1.007276466879) * self.charge

    def getFragMethod(self):
        return self.fragMethod

    def getDetector(self):
        return self.detector

    def printf(self):
        print ("Title, RT, PEPMASS, PEPINT, CHARGE")
        print (self.title, self.RT, self.pepmz, self.pepint, self.charge)

    def to_mgf(self):
        # need dummy values in case no peak charges are in the data
        if len(self.peakcharge) == 0:
            self.peakcharge = [""]*self.peaks.shape[0]
        mgf_str = """
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (self.title, self.RT, self.pepmz, self.pepint, self.charge, "\r\n".join(["%s %s %s" % (i[0], i[1], j, ) for i,j in zip(self.peaks, self.peakcharge)]))
        return mgf_str

    
#==============================================================================
# File Reader
#==============================================================================
class MGF_Reader():
    """A MGF_Reader is associated with a FASTA file or an open connection
    to a file-like object with content in FASTA format.
    It can generate an iterator over the sequences.

    Usage:
    --------------------------
    >>reader = MGF_Reader() \r\n
    >>reader.load(infile) \r\n
    >>#do something \r\n
    >>reader.store(outfile, outspectra) \r\n
    """
    def __init__(self, verbose=False):
        self.verbose = verbose
        
    def load(self, infile, getpeakcharge=False):
        """
        Function to set the input file for the MGF file.

        Parameters:
        -----------------------------
        infile: str,
                file location



        """
        self.infile = infile
        self.peakcharge = getpeakcharge

    def __iter__(self):
        mgf_file = open(self.infile, "r")
        found_ions = False
        for line in mgf_file:
            if len(line.strip()) == 0:
                continue
            if line.startswith("BEGIN IONS"):
               # init arrays for peaks
                found_ions = True
                mass = []
                intensity = []
                peakcharge = []
            elif line.startswith("TITLE"):
                title = re.search("TITLE=(.*)", line).groups()[0]

            elif line.startswith("RTINSECONDS"):
                RT = float(re.search("RTINSECONDS=(.*)", line).groups()[0])

            elif line.startswith("PEPMASS"):
                precursor = re.search("PEPMASS=(.*)", line).groups()[0].split()
                pep_mass = float(precursor[0])
                try:
                    pep_int = float(precursor[1])
                except:
                    pep_int = -1.0

            elif line.startswith("CHARGE"):
                charge = float(re.search("CHARGE=(\d)", line).groups()[0])

            elif "=" in line:
                if self.verbose:
                    print ("unhandled paramter: %s" % (line))

            elif line.startswith("END IONS"):
                if sys.version_info.major < 3:
                    ms = MS2_spectrum(title, RT, pep_mass, pep_int, charge, np.array(zip(mass, intensity)), peakcharge)
                else:
                    ms = MS2_spectrum(title, RT, pep_mass, pep_int, charge, np.array((mass, intensity)), peakcharge)

                yield ms
            else:
                if found_ions is True:
                    peak = line.split()
                    mass.append(float(peak[0]))
                    intensity.append(float(peak[1]))
                    if self.peakcharge:
                        if len(peak) >2:
                            peakcharge.append(peak[2])

    def store(self, out_file, ms_list, write_charge=False):
        """
        Function to store data as MGF file.

        Parameters:
        ------------------------
        out_file: str,
                  output file location
        ms_list: list or iterable
                 container of ms2 objects that are writte to the file.
        """
        out_mgf = open(out_file, "w")

        for ms in ms_list:
            if write_charge:
                mgf_str = """
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
""" % (ms.title, ms.RT, ms.pepmass, ms.pepint, ms.charge, "\r\n".join(["%s %s %s" % (i[0], i[1], j, ) for i,j in zip(ms.peaks, ms.peakcharge)]))
            else:
                mgf_str = """
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
""" % (ms.title, ms.RT, ms.pepmass, ms.pepint, ms.charge, "\r\n".join(["%s %s" % (i, j) for i,j in ms.peaks]))

            out_mgf.write(mgf_str)
        out_mgf.close()



class APL_Reader():
    """A APL_Reader is associated with a apl (MaxQuant) file or an open connection
    to a file-like object with content in (MaxQuant) format.
    It can generate an iterator over the sequences.

    Usage:
    --------------------------
    >>reader = MGF_Reader() \r\n
    >>reader.load(infile) \r\n
    >>#do something \r\n
    >>reader.store(outfile, outspectra) \r\n
    """
    def load(self, infile):
        """
        Function to set the input file for the MGF file.

        Parameters:
        -----------------------------
        infile: str,
                file location

        Usage:
        --------------------------
        >>reader = APL_Reader() \r\n
        >>reader.load(infile) \r\n
        >>#do something \r\n
        >>reader.store(outfile, outspectra) \r\n

        """
        self.infile = infile

    def __iter__(self):
        mgf_file = open(self.infile, "r")
        found_ions = False
        for line in mgf_file:
            if len(line.strip()) == 0:
                continue
            if line.startswith("peaklist start"):
               # init arrays for peaks
                found_ions = True
                mass = []
                intensity = []
                peakcharge = []
            elif line.startswith("header"):
                title = re.search("header=(.*)", line).groups()[0]
            elif line.startswith("charge"):
                charge = float(re.search("charge=(\d)", line).groups()[0])

            elif line.startswith("mz"):
                # actually izs MZ!
                pep_mass = float(re.search("mz=(.*)", line).groups()[0])

            elif line.startswith("fragmentation"):
                    continue

            elif "=" in line:
                print ("unhandled paramter: %s" % (line))

            elif line.startswith("peaklist end"):
                ms = MS2_spectrum(title, -1, pep_mass, -1, charge, np.array(zip(mass, intensity)), peakcharge)
                yield ms
            else:
                if found_ions is True:
                    peak = line.split("\t")
                    if len(peak) == 1:
                        peak = line.split()
                    mass.append(float(peak[0]))
                    intensity.append(float(peak[1]))
                    if len(peak) >2:
                        peakcharge.append(float(peak[2]))

    def store(self, out_file, ms_list):
        """
        Function to store data as MGF file.

        Parameters:
        ------------------------
        out_file: str,
                  output file location
        ms_list: list or iterable
                 container of ms2 objects that are writte to the file.
        """
        out_mgf = open(out_file, "w")
        for ms in ms_list:

            mgf_str = """
peaklist start
header=%s
mz=%s
charge=%s
%s
peaklist end
""" % (ms.title, ms.pepmass, ms.charge, "\r\n".join(["%s %s" % (i, j) for i,j in ms.peaks]))
            out_mgf.write(mgf_str)
        out_mgf.close()
