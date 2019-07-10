#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 17:50:25 2018

@author: sven
"""

import pyopenms as oms
#from PySpecViewer import TheoreticalSpectrumGenerator as TSG
import numpy as np
import re
import matplotlib.pyplot as plt

PROTON = 1.00727646658

def find_parent_MS1_scan(exp, scan):
    """
    Gets the scan index of the preceding MS1 scan.
    """
    counter = 1
    while True:
        if exp[scan-counter].getMSLevel() == 1:
            return(scan-counter)
        else:
            counter += 1
        
def get_error(mz_calc, mz_exp):
    """
    Computes the absolute error in ppm
    """
    return((np.abs((mz_calc - mz_exp)) / mz_exp)*10**6)

def get_scan(nativeid):
    """
    Retrieves the scan number from the native header id.
    
    Parameters:
    ------------------
    nativeid: str,
                header of the MGF (mzML)
    """
    return(int(re.search("scan=(\d+)", nativeid).groups()[0]))
    
def extend_mass_mz(exp, MS1scan, mz, seed_scan, scans_masstrace, 
                   peaks_masstrace, ppm_trace, add=True, RTdiff_max=90):
    """
    Given an m/z; tries to find the same m/z along the RT dimension.
    
    Parameters:
    --------------
    exp: oms exp,
            experiment, parsed mzML file
            
    MS1scan: int,
            scan number
    mz: float,
        expected m/z of the precursor
        
    seed_scan: int,
               the MS2 scan
               
    scans_masstrace: ar-like,
                store results here
        
    peaks_masstrace: ar-like,
                    store pekas here
    ppm_trace: float,
                ppm error to tolerate
    add: bool,
            if True go into + RT direction, if False go into - RT direction
    
    RTdiff: float,
            determines how far to look for the given m/z precursor
    """
    if add==True:
        constant = 1
    else:
        constant = -1
        
    RTdiff = 0
    scan_it = -1
    maxcount = 200
    currentcount = 0
    while RTdiff <= RTdiff_max:
        currentcount += 1
        if currentcount >= maxcount or (MS1scan + scan_it <= 1):
            # print ("max nscan diff reached") # TODO: check if really sensible (data I checked seemed ok for a bit over 200)
            break
        scan_it += constant
        temp_scan = exp[MS1scan + scan_it]
        RTdiff = np.abs(seed_scan.getRT() - temp_scan.getRT())
        if temp_scan.getMSLevel() == 1:
            candidate = temp_scan.findNearest(mz)
            if get_error(temp_scan[candidate].getMZ(), mz) <= ppm_trace:
                scans_masstrace.append(get_scan(str(temp_scan.getNativeID())))
                peaks_masstrace.append(temp_scan[candidate])
    
def extend_isotope_mz():
    """
    Given a m/z peak; tries to find the isotope clusters that belong to that 
    peak.
    """
    #TODO
    
def extract_mass_trace(exp, MS1scan, mz, charge, ppm, ppm_trace, RTdiff=90):
    """
    Extracts a summed ion chromatogram for a given m/z and the isotope peaks
    
    
    Parameters:
    
    """
    #set the seed data, e.g. which scan was the MS2, which peak is the closest
    # to the given precursor (mz)
    seed_scan = exp[MS1scan]
    seed_peak = seed_scan.findNearest(mz)
    diff = PROTON / charge
    
    #get the peaks of the seed scan and arrange into matrix
    peaks = seed_scan.get_peaks()
    peaks = np.column_stack(peaks)
    peak_clust = 1
    seeds = [peaks[seed_peak][0]]
    
    #collect the scans, and peaks in the mass trace
    scans_masstrace = [MS1scan]
    peaks_masstrace = []
    
    #first extend the precursor mz in +RT direction
    extend_mass_mz(exp, MS1scan, mz, seed_scan, scans_masstrace, 
                   peaks_masstrace, ppm_trace, RTdiff_max=RTdiff, add=True)
    
    #then extend the precursor mz in -RT direction
    extend_mass_mz(exp, MS1scan-1, mz, seed_scan, scans_masstrace, 
                   peaks_masstrace, ppm_trace, RTdiff_max=RTdiff, add=False)
            
    peaks_mt = np.matrix([(i.getMZ(), i.getIntensity()) for i in peaks_masstrace])
    return(peaks_mt, scans_masstrace)
    #TODO extend isotope
#    
#    plt.plot(peaks_mt[:,1], '--o')
#    
#    #find the most abundant ion and use as seed for isotope extension
#    best_isotope_seed = np.argmax(peaks_mt[:,1])
#    isotope_seed_scan = exp[scans_masstrace[best_isotope_seed]]
#    istp_seed_peaks = np.column_stack(isotope_seed_scan.get_peaks())
#    istp_seed_peak = isotope_seed_scan.findNearest(mz)
#    istp_seeds = []
#    for i in range(0, 5):
#        expected = istp_seed_peaks[istp_seed_peak + i-1][0] \
#                    + (PROTON / charge) * peak_clust
#        next_peak = istp_seed_peaks[istp_seed_peak + i][0]
#        if get_error(next_peak, expected) <= ppm:
#            istp_seeds.append(next_peak)
#            
#            
#    plt.bar(peaks[275:280][:,0], peaks[275:280][:,1], width=0.03)
#    
#    plt.bar(peaks[:,0], peaks[:,1])

# =============================================================================
# test
# =============================================================================
#example

if __name__ == '__main__':
    mzml_file = "data/MS1_data/E130808_16_Quty1_AB_DE_220_V127_H_6_666ul.mzML"
    #scan_dictionary = read_mzmls(path="data/MS1_data/E130808_16_Quty1_AB_DE_220_V127_H_6_666ul")

    mzml = oms.MzMLFile()
    exp = oms.MSExperiment()
    mzml.load(mzml_file, exp)


    #init dictionary
    basename = "Test"
    scan_dictionary = {basename:{}}
    #iteate over spectra
    for spectrum in exp:
        if spectrum.getMSLevel() == 2:
            scan = int(re.search("scan=(\d+)", str(spectrum.getNativeID())).groups()[0])
            scan_dictionary[basename][scan] = spectrum

    #set example
    #we need:
    # exp, MS1scan, mz, seed_scan, scans_masstrace, peaks_masstrace, ppm_trace, add=True, RTdiff=120
    MS2scan = 30018
    MS1scan = find_parent_MS1_scan(exp, MS2scan)
    precursor = exp[MS2scan].getPrecursors()[0]

    mz = precursor.getMZ()
    charge = precursor.getCharge()
    ppm = 10
    ppm_trace = 5
    RT = exp[MS2scan].getRT()
    print (mz, charge, MS1scan, MS2scan, RT, RT/60.)

    #extract the mass trace!
    mz_trace, scans_trace = extract_mass_trace(exp, MS1scan, mz, charge, ppm, ppm_trace, 60)
    #use the peak with the highest intensity as seed
    best_isotope_seed = np.argmax(mz_trace[:,1])
    #get the spectrum where the peak had the highest intensity
    best_seed_spectrum = scans_trace[best_isotope_seed]
    plt.plot(mz_trace[:,1], '-')
    plt.xlabel("nscans")
    plt.ylabel("Intensity")

#gasphase stuff       
#newscan = scan-1
#pep1_mass = get_peptide_mass(TSG, pepseq1)
#pep2_mass = get_peptide_mass(TSG, pepseq2)
#pep1_mass_acl = pep1_mass + cl_mass
#pep2_mass_acl = pep2_mass + cl_mass
#cl_mz = raw_scan.getPrecursors()[0].getMZ()
#cl_charge = raw_scan.getPrecursors()[0].getCharge()
#RT = raw_scan.getRT()


