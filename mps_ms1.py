import os
import pyopenms as oms
import numpy as np
from ProteoFileReader import MGF_Reader
import mass_trace
from joblib import Parallel, delayed
import re


def add_relaxation_mgf(mgf, mps, outfile, create_comparison=False):
    # mgf_file, outdir, differences = argls[0], argls[1], argls[2]
    mass_diff = 1.00335483
    # if '.mgf' in args['file']:
    filename = os.path.split(mgf)[1]
    # read mgf
    spectra = MGF_Reader()
    spectra.load(mgf)
    out_writer = open(outfile, "w")
    for spectrum in spectra:
        # calculate mass (neglect proton bec. later on difference used)
        regex_match = re.search('(scan=)[0-9]*', spectrum.getTitle())
        if regex_match is not None:
            scan = int(regex_match.group(0).split('scan=')[1])
        else:
            scan = int(spectrum.getTitle().split('.')[-2])
        # scan = int(spectrum.getTitle().split('.')[-2])
        # try:
        differences = [0, -1, -2, -3, -4]
        if not create_comparison:
            row = mps[mps[:, 0] == scan, 1:]
            if len(row) == 1:
                differences = [-i for i in range(len(row[0])) if row[0][i] == 1]
            elif len(row) > 1:
                raise ValueError('multiple matches to scan %s' % scan)
            else:
                print 'scan %s not found' % scan

        mass = spectrum.getPrecursorMZ() * spectrum.charge          # ToDo: correct? talk to SL -protonMass missing!
        spectra_add_mip = [str((mass + x * mass_diff) / spectrum.charge) for x in differences if x != 0]
        if 0 in differences:
            prec_mz = spectrum.getPrecursorMZ()
        else:
            prec_mz = spectra_add_mip[0]
            spectra_add_mip = spectra_add_mip[1:]
        # except KeyError:
        #     differences = [-2, -1, 0]

        stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
ADDITIONALMZ={}
{}
END IONS     """.format(spectrum.getTitle(), prec_mz,
                        spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                        int(spectrum.charge), spectrum.getRT(),
                        ';'.join(spectra_add_mip),
                "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks]))
        out_writer.write(stavrox_mgf)


def read_mzml(infile):
    # init variables
    mzml_file = oms.MzMLFile()
    exp = oms.MSExperiment()

    # load spectra into exp
    mzml_file.load(infile, exp)
    return (exp)


def get_error(mz1, mz2, charge=None, ppm=True):
    if ppm:
        return abs(mz1 - mz2) / mz2 * 1e6
    else:
        if charge is not None:
            return (mz1 - mz2) * charge


def return_mps_range(intensity, mps_max):
    if False: #intensity >= 10e6:
        return [-1]
    else:
        return list(range(-1, mps_max - 1, -1))


def ms1_peaks(exp, tolerance=6, mps_range=[-1, -2, -3, -4]): # mps_range=[-1, -2, -3, -4]
    # loop through spectra and make count if precursor in MS2 spectrum
    peaks_found = []
    n_ms2 = 0
    nspectra = exp.size()
    for i, spectrum in enumerate(exp):

        if spectrum.getMSLevel() == 1:
            continue
        n_ms2 += 1
        if i % 10000 == 0:
            print "{}/{} done..".format(i, nspectra)

        # if iidone not in matched_scans:
        #     continue
        precursor = spectrum.getPrecursors()[0]
        prec_mz, prec_charge = precursor.getMZ(), precursor.getCharge()
        MS1scan = mass_trace.find_parent_MS1_scan(exp, i)
        ppm_pseudo = 10 # taken from Svens script, apparently not used in function
        mz_trace, scans_trace = mass_trace.extract_mass_trace(exp, MS1scan, prec_mz, prec_charge, ppm_pseudo, tolerance, 10) # test with 20 to see difference
        if len(scans_trace) == 1:
            tmp_ms1 = exp[MS1scan] # TODO find out why mass_trace does not find anything
        else:
            # try:
            best_isotope_seed = np.argmax(mz_trace[:, 1])
            # except:
            #     pass
            best_seed_spectrum = scans_trace[best_isotope_seed]

            # ms1_prev = 0
            # for j in range(i, 0, -1):
            tmp_ms1 = exp[best_seed_spectrum]

        res = tmp_ms1[tmp_ms1.findNearest(prec_mz)]
        if abs(res.getMZ() - prec_mz) / prec_mz <= tolerance:
            prec_int = res.getIntensity()
            # mps_range = return_mps_range(prec_int, mps_max)
            # if prec_int >= 3e6:
            #     peaks_found.append([i + 1, True] + [False] * len(mps_range))
            #     continue
            # else:
            #     peaks_found.append([i + 1, True] + [True] * len(mps_range))
            #     continue
            theo_mip = np.array([prec_mz + (mip_i * 1.00335483) / prec_charge for mip_i in mps_range])
            mip_nearest = [tmp_ms1.findNearest(x) for x in theo_mip]
            error = np.array([get_error(tmp_ms1[mea].getMZ(), expi, ppm=True) for expi, mea in
                              zip(theo_mip, mip_nearest)])
            range_found = [True if x <= tolerance else False for x in error]
            if sum(range_found) == 0:
                # TODO: try if sensible to not mps search these
                if len(mps_range) > 2:
                    peaks_found.append(
                        # [i + 1, True] + [True] * len(mps_range)
                        # [i + 1, True, True, True] + [False] * (len(mps_range) - 2)
                        # [i + 1, True, True, True, True, False]
                        [i + 1, True] + [True] * (len(mps_range) - 1) + [False]
                    )
                else:
                    peaks_found.append(
                        [i + 1, True] + [True] * len(mps_range)
                    )
                continue
            else:
                # check for continous peaks except -1 peak
                # TODO: allow gap?
                found = False
                lightest_peak = 1
                for i_mip in range(len(range_found), 1, -1):
                    if range_found[i_mip - 1] & (i_mip > lightest_peak):
                        lightest_peak = i_mip
                    if sum(range_found[:i_mip]) == len(range_found[:i_mip]):
                        # if i_mip == len(range_found):
                        #     sel = [False] * (i_mip - 2) + range_found[i_mip - 2:]
                        # else:
                        #     sel = [False] * (i_mip - 2) + [True] * 3 + [False] * (len(range_found) - 1 - i_mip) # 2
                        # takes lightest 2 continous + existing lighter peaks, excludes heaviar
                        sel = [False] * (i_mip - 2) + range_found[i_mip - 2:] # 2
                        peaks_found.append([i + 1, False] + sel)
                        found = True
                        break
                # if no continous found take all
                if not found:
                    if not lightest_peak == len(range_found):
                        peaks_found.append(
                            # [i + 1, True] + [True] * len(mps_range)
                            # [i + 1, True] + [True] * (lightest_peak + 1) + [False] * (len(range_found) - lightest_peak - 1)
                            [i + 1, True] + [True] * (lightest_peak) + [False] * (len(range_found) - lightest_peak)
                        )
                    else:
                        peaks_found.append([i + 1, True] + [True] * len(mps_range))
                    continue

        else:
            print 'Precursor not found'
            continue

    return np.array(peaks_found)


def main(mzmlfile, exp_id, setting, infoout_dir, mgf_in_dir, mgf_out_dir):
    exp = read_mzml(mzmlfile)
    # exp_id = mzml_file[:10]

    mps_df = ms1_peaks(exp)
    np.savetxt(infoout_dir + '/%s_%s.csv' % (setting, exp_id), mps_df, delimiter=',')
    corresponding_mgf = [x for x in os.listdir(mgf_in_dir) if exp_id in x][0]
    add_relaxation_mgf(mgf=mgf_in_dir + corresponding_mgf, mps=mps_df,
                       outfile=mgf_out_dir + '/%s_' % setting + corresponding_mgf)

if __name__ == '__main__':
    isotope_diff = 1.00335483
    # mzml_dir = 'D:/user/Swantje/data/PC/mzML/'
    # chaet_dir = 'fr7-10'
    mzml_dir = 'D:/user/Swantje/dsso/myco_prepro/1c_remaining_frac/test_set/mzml/'
    mgf_filtered_dir = '//130.149.167.198/rappsilbergroup/users/lswantje/dsso_opt/searches/mycoplasma_opt/data/filtered/'
    # mgf_filtered_dir = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/%s/All_prepro_peakfiles/mscon_PF_20/' % chaet_dir
    setting_name = 'mpsreductionfilt'
    mgf_out = '//130.149.167.198/rappsilbergroup/users/lswantje/dsso_opt/searches/mycoplasma_opt/data/' + setting_name
    info_out = '//130.149.167.198/rappsilbergroup/users/lswantje/dsso_opt/searches/mycoplasma_opt/' + setting_name
    # mgf_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/%s/All_prepro_peakfiles/' % chaet_dir + setting_name
    # info_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/%s/relaxation_tbls/' % chaet_dir + setting_name

    if not os.path.exists(mgf_out):
        os.makedirs(mgf_out)
    if not os.path.exists(info_out):
        os.makedirs(info_out)

    mzmls_in = [x for x in os.listdir(mzml_dir) if '.mzML' in x]
    # Parallel(n_jobs=4)(delayed(main)(mzml_dir + x, x[:10], setting_name, info_out, mgf_filtered_dir,
    #                                  mgf_out) for x in mzmls_in)
    for x in mzmls_in:
        main(mzml_dir + x, x[:10], setting_name, info_out, mgf_filtered_dir, mgf_out)
    for mzml_file in [x for x in os.listdir(mzml_dir) if '.mzML' in x]:
        exp = read_mzml(mzml_dir + mzml_file)
        exp_id = mzml_file[:10]

        mps_df = ms1_peaks(exp)
        np.savetxt(info_out + '/%s_%s.csv' % (setting_name, mzml_file), mps_df, delimiter=',')
        corresponding_mgf = [x for x in os.listdir(mgf_filtered_dir) if exp_id in x][0]
        add_relaxation_mgf(mgf=mgf_filtered_dir + corresponding_mgf, mps=mps_df,
                           outfile=mgf_out + '/%s_' % setting_name + corresponding_mgf)

    # for mgf_file in os.listdir(mgf_filtered_dir):
    #     add_relaxation_mgf(mgf=mgf_filtered_dir + mgf_file,
    #                        outfile='D:/user/Swantje/projects/pipeline_prepro_xi_fdr/lars_PC_4frag_BS3_Lumos/All_prepro_peakfiles/' + setting_name + '_' + mgf_file)
