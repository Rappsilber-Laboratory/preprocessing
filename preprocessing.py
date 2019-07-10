import os
import subprocess
from multiprocessing import Pool
import sys
import re
import getopt
from pyteomics import mzml
from functools import partial
import ProteoFileReader
import mass_recal_ms2
import mass_recal


def read_cmdline():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['input=', 'config=', 'outpath=', 'db=', 'xiconf=', 'shiftcsv=',
                                                      'skip_recal=', 'skip_ms2_recal='])
    except getopt.GetoptError:
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file> '
              '--db <path to database to search for recalibration>'
              '--xiconf <path to xi config to use for recalibration>',
              '--shiftcsv <path to csv with fixed shifts>',
              '--skip_recal <boolean>',
              '--skip_ms2_recal <boolean>')
        sys.exit()
    recal = True
    ms2recal = True
    recal_conf = {}
    for opt, arg in opts:
        if opt == '--input':
            input_arg = arg
        elif opt == '--outpath':
            outdir = arg
        elif opt == '--config':
            config = arg
        elif opt == '--db':
            recal_conf['db'] = arg
        elif opt == '--xiconf':
            recal_conf['xiconf'] = arg
        elif opt == '--shiftcsv':
            recal_conf['shift_csv'] = arg
        elif opt == '--skip_recal':
            recal = False
        elif opt == '--skip_ms2_recal':
            ms2recal = False

    if 'input_arg' not in locals() or 'config' not in locals():
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file> '
              '--db <path to database to search for recalibration>'
              '--xiconf <path to xi config to use for recalibration>',
              '--shiftcsv <path to csv with fixed shifts> --skip_recal <boolean>')
        sys.exit()
    # if no outdir defined use location of input
    if 'outdir' not in locals() and os.path.isdir(input_arg):
        outdir = os.path.join(input_arg, 'processed')
    elif 'outdir' not in locals() and not os.path.isdir(input_arg):
        outdir = os.path.join(os.path.split(input_arg)[0], 'processed')
    if 'shift_csv' not in recal_conf:
        recal_conf['shift_csv'] = None
    if ('db' not in recal_conf or 'xiconf' not in recal_conf) and recal:
        print('Recalibration enabled but parameters missing! Set --db and --xiconf.'
              'preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file> '
              '--db <path to database to search for recalibration>'
              '--xiconf <path to xi config to use for recalibration>',
              '--shiftcsv <path to csv with fixed shifts> --skip_recal <boolean>')
        sys.exit()

    return input_arg, outdir, config, recal_conf, recal, ms2recal


def mzml_to_MS2_spectra(mzml_file, detector_filter="all"):
    """
    function to split a mzML file into a list of MS2_Spectra objects (can be written to mgf format)
    with fragmentation method and detector type

    Parameters:
    -----------------------------------------
    mzml_file: str,
            path to mzML file
    detector_filter: filter scans by detector type ('all', 'FT', 'IT')

    Return: list(MS2_spectrum)

    """

    mzml_reader = mzml.read(mzml_file)
    sorted_ms2_spectra = []
    unknown_frag_method_count = 0

    n = 0
    for spectrum in mzml_reader:
        if spectrum['ms level'] == 2:
            n += 1
            filter_str = spectrum['scanList']['scan'][0]['filter string']
            try:
                detector_str = re.search("^(FT|IT)", filter_str).groups()[0]
                frag_groups = re.findall("@([A-z]+)([0-9.]+)", filter_str)
            except AttributeError:
                raise Exception("filter string parse error: %s" % filter_str)

            if not detector_filter == "all":
                if not detector_filter == detector_str:
                    continue

            title = os.path.split(mzml_file)[1].split('.mzML')[0] + " " + spectrum['id']
            rt = spectrum['scanList']['scan'][0]['scan start time'] * 60
            precursor = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            pre_mz = precursor['selected ion m/z']
            try:
                pre_int = precursor['peak intensity']
            except KeyError:
                pre_int = 0
            pre_z = precursor['charge state']
            peaks = zip(spectrum['m/z array'], spectrum['intensity array'])

            frag_methods = [f[0] for f in frag_groups]

            if "etd" in frag_methods:
                if "cid" in frag_methods:
                    frag_method = "ETciD"
                elif "hcd" in frag_methods:
                    frag_method = "EThcD"
                else:
                    frag_method = "ETD"
            elif "cid" in frag_methods:
                frag_method = "CID"
            elif "hcd" in frag_methods:
                frag_method = "HCD"
            else:
                frag_method = 'unknown'
                unknown_frag_method_count += 1

            ms2class_spectrum = ProteoFileReader.MS2_spectrum(
                title,
                rt,
                pre_mz,
                pre_int,
                pre_z,
                peaks,
                detector=detector_str,
                fragmethod=frag_method
            )

            sorted_ms2_spectra.append(ms2class_spectrum)

    if unknown_frag_method_count > 0:
        raise Warning("The fragmentation method of %i spectra could not be identified" % unknown_frag_method_count)

    return sorted_ms2_spectra


def generate_cihcd_spectra(mzml_file):
    """

    """

    mzml_reader = mzml.read(mzml_file)
    cihcd_spectra = []

    n = 0
    for spectrum in mzml_reader:
        if spectrum['ms level'] == 3:
            n += 1
            filter_str = spectrum['scanList']['scan'][0]['filter string']
            try:
                detector_str = re.search("^(FT|IT)", filter_str).groups()[0]
                frag_groups = re.findall("@([A-z]+)([0-9.]+)", filter_str)
                precursor_mz_groups = re.findall("([0-9.]+)@", filter_str)
            except AttributeError:
                raise Exception("filter string parse error: %s" % filter_str)
            try:
                ms2_id = spectrum['precursorList']['precursor'][0]['spectrumRef']
            except KeyError:
                ms2_id = ''  # TODO why Key Error
            title = os.path.split(mzml_file)[1].split('.mzML')[0] + " " + spectrum['id'] + " ms2_scanId=" + ms2_id
            rt = spectrum['scanList']['scan'][0]['scan start time'] * 60

            pre_mz = precursor_mz_groups[0]     # take ms2 precursor as precursor
            pre_int = -1
            pre_z = -1
            peaks = zip(spectrum['m/z array'], spectrum['intensity array'])

            ms2class_spectrum = ProteoFileReader.MS2_spectrum(title, rt, pre_mz, pre_int, pre_z, peaks)

            cihcd_spectra.append(ms2class_spectrum)

    return cihcd_spectra


def mscon_cmd(filepath, outdir, settings, mgf):
    filename = os.path.split(filepath)[1]

    if (filename[:filename.rfind('.')] in [x[:x.rfind('.')] for x in os.listdir(outdir)]) or (os.path.isdir(filepath)):
        print('File ' + filename + ' already existing in output directory, not done again.')
        return []

    filter_formatted = []
    for i in range(len(settings)):
        filter_formatted.append('--filter')
        filter_formatted.append(settings[i])

    if mgf:
        cmd_list = [filepath, '--mgf', '-o', outdir] + filter_formatted
    else:
        cmd_list = [filepath, '-o', outdir] + filter_formatted
    return cmd_list


def process_file(filepath, outdir, mscon_settings, split_acq, detector_filter, mscon_exe, cihcd_ms3=False):  #TODO implement option further up
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    conv_cmds = mscon_cmd(filepath=filepath, outdir=outdir, settings=mscon_settings, mgf=not split_acq)

    if len(conv_cmds) > 0:
        msconvert = subprocess.Popen([mscon_exe] + conv_cmds)
        msconvert.communicate()

    filename = os.path.split(filepath)[1]
    mzml_file = os.path.join(outdir, filename[:filename.rfind('.')] + '.mzML')

    if cihcd_ms3:
        cihcd_spectra = generate_cihcd_spectra(mzml_file)
        ProteoFileReader.write_mgf(spectra=cihcd_spectra, outfile=os.path.join(outdir, 'CIhcD_ms3_' + filename[:filename.rfind('.')] + '.mgf'))

    if split_acq:
        split_spectra = mzml_to_MS2_spectra(mzml_file, detector_filter)

        ProteoFileReader.write_mgf(
            spectra=split_spectra,
            outfile=os.path.join(outdir, filename[:filename.rfind('.')]+'.mgf')
        )


if __name__ == '__main__':
    # read cmdline arguments / get default values
    input_arg, outdir, config_path, recal_conf, recal, ms2recal = read_cmdline()
    try:
        execfile(config_path)
    except NameError:
        exec(open(config_path).read())

    # get files in directory
    if os.path.isdir(input_arg):
        full_paths = [os.path.join(input_arg, rawfile) for rawfile in os.listdir(input_arg) if rawfile[-4:] == '.raw']
        full_paths = [x for x in full_paths if not os.path.isdir(x)]
    # if single file given reformat to list
    # TODO allow txt file with filepath
    else:
        full_paths = [input_arg]

    print("""file input:
{}
""".format('\n'.join(full_paths)))

    # start msconvert for conversion and peak filtering
    pool = Pool(processes=nthr)
    pool.map(partial(process_file, outdir=outdir, mscon_settings=mscon_settings, split_acq=split_acq,
                     detector_filter=detector_filter, mscon_exe=msconvert_exe), full_paths)
    pool.close()
    pool.join()

    mgf_file_list = [os.path.join(outdir, x) for x in os.listdir(outdir) if '.mgf' in x]
    if recal:

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # TODO change to parallel with manual input of error
        for inputfile in mgf_file_list:
            if 'ms3' in os.path.split(inputfile)[1]:
                continue
            if ms2recal:
                mass_recal_ms2.main(fasta=recal_conf['db'], xi_cnf=recal_conf['xiconf'], outpath=outdir,
                                    mgf=inputfile, threads=str(nthr),
                                    val_input=recal_conf['shift_csv']
                                    )
            else:
                mass_recal.main(fasta=recal_conf['db'], xi_cnf=recal_conf['xiconf'], outpath=outdir,
                                mgf=inputfile, threads=str(nthr),
                                val_input=recal_conf['shift_csv']
                                )

        mgf_file_list = [os.path.join(os.path.split(x)[0], 'recal_' + os.path.split(x)[1]) for x in mgf_file_list]

    if split_acq:
        for mgf_file in mgf_file_list:
            ProteoFileReader.split_mgf_methods(mgf_file)
