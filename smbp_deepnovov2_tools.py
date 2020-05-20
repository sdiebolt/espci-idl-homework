# -*- coding: utf-8 -*-
"""Create DeepNovoV2 input files from MGF spectra and Mascot annotation files.

This module contains tools to annotate MGF spectra using Mascot XML results and
create the features.csv file needed by DeepNovoV2.
"""

import pandas as pd
import pyteomics.mgf as mgf
import re
from typing import List
import xml.etree.ElementTree as ET


def extract_features(mgf_input: str, mascot_input: str = None) -> pd.DataFrame:
    """Extract features from an MGF file to create DeepNovoV2 features.csv.

    If a Mascot XML results file is provided, sequences will be added to the
    features and non-identified spectra will be discarded.

    Args:
        mgf_input (str): path to the MGF file.
        mascot_input (str): path to the Mascot XML file.

    Returns:
        pandas.DataFrame: the features DataFrame.

    """
    # Extract mascot sequences if a mascot XML file was given.
    if mascot_input:
        mascot_seq = extract_mascot_sequences(mascot_input)

    with mgf.read(mgf_input, read_charges=False) as reader:
        rows_list = list()
        for spectrum in reader:
            # Retrieve spectrum parameters
            params = spectrum['params']
            scans = params['scans']

            # In Mascot MGF files, the PEPMASS line contains the precursor m/z,
            # intensity and charge: only m/z is selected.
            mz = (params['pepmass'][0] if len(params['pepmass']) > 1 else
                  params['pepmass'])
            # Pyteomics returns charges as a list of integers with subtype
            # pyteomics.auxiliary.structures.Charge. The first charge is
            # selected and converted to integer.
            charge = int(params['charge'][0])

            rt_mean = params['rtinseconds']

            # Get sequences for the current spectrum if a mascot XML file was
            # given. If no XML file was provided, add an empty sequence.
            if mascot_input:
                sequences = mascot_seq.loc[
                    mascot_seq.title == params['title'], 'sequence'
                ].values
            else:
                sequences = ['']

            # Add one feature per sequence found.
            for seq in sequences:
                # spec_group_id and scans are identical (DeepNovoV2 uses these
                # columns in case of merged spectra).
                # profile and feature area are set to default values set by the
                # authors.
                rows_list.append({
                    'spec_group_id':  scans,
                    'scans': scans,
                    'm/z': mz,
                    'z': charge,
                    'rt_mean': rt_mean,
                    'seq': seq,
                    'profile': '0.0:1.0',
                    'feature area': '1.0'
                })

    # Return DataFrame of features.
    return pd.DataFrame(rows_list).drop_duplicates()


def annotate_mgf(mgf_input: str, mascot_input: str, mgf_output: str):
    """Annotate MGF file using Mascot XML results.

    annotate_mgf will annotate the MGF file using peptide sequences found in
    the Mascot XML results and write the resulting MGF file to mgf_output.

    Args:
        mgf_input (str): path to the MGF input file.
        mascot_intput (str): path to the Mascot XML results.
        mgf_output (str): path to the MGF output file.

    """
    # Retrieve mascot sequences.
    mascot_seq = extract_mascot_sequences(mascot_input)

    with mgf.read(mgf_input, read_charges=False) as reader:
        for spectrum in reader:
            sequences = mascot_seq.loc[
                mascot_seq.title == spectrum['params']['title'], 'sequence'
            ].values
            # If multiple sequences are associated to a single spectrum, the
            # latter will be duplicated for each sequence.
            for seq in sequences:
                spectrum['params']['seq'] = seq
                mgf.write((spectrum,), mgf_output)


def merge_mgf(input_files: List[str], output_file: str):
    """Merge a list of MGF files, re-numbering the scan IDs.

    Args:
        input_files (list): a list of paths to MGF files.
        output_file (str): path to the merged MGF file.

    """
    with open(output_file, mode="w") as output_handle:
        for index, input_file in enumerate(input_files):
            with open(input_file, mode="r") as input_handle:
                line = input_handle.readline()
                while line:
                    # If the current line is a scan ID, renumber it. Otherwise,
                    # write the line as is.
                    if "SCANS=" in line:
                        scan = re.split('=|\n', line)[1]
                        output_handle.write("SCANS=F{0}:{1}\n"
                                            .format(index, scan))
                    else:
                        output_handle.write(line)
                    line = input_handle.readline()


def merge_features(input_files: List[str], output_file: str):
    """Merge a list of feature files, re-numbering the scan IDs.

    Args:
        input_files (list): a list of paths to features files.
        output_file (str): path to the merged MGF file.

    """
    features_list = list()
    for index, f in enumerate(input_files):
        features = pd.read_csv(f)
        # Renumber the spectrum group and scan IDs.
        features.spec_group_id = ('F' + str(index) + ':' +
                                  features.spec_group_id.astype(str))
        features.scans = ('F' + str(index) + ':' +
                          features.scans.astype(str))
        features_list.append(features)
    pd.concat(features_list).to_csv(output_file, index=False)


def extract_mascot_sequences(mascot_xml: str):
    """Parse Mascot XML results and retrieve sequences with modifications.

    Only first rank peptides will be kept.

    Args:
        mascot_xml (str): path to the Mascot XML result file.

    Returns:
        pandas.DataFrame: DataFrame containing spectrum titles and
            their sequences.

    """
    ns = {'m':
          'http://www.matrixscience.com/xmlns/schema/mascot_search_results_2'}
    mascot_root = ET.parse(mascot_xml).getroot()

    # Extract variable modifications identifiers.
    modifications = mascot_root.findall('m:variable_mods/m:modification', ns)
    mods = dict()
    for m in modifications:
        mods.update(
            {m.get('identifier'):
             re.sub(' (.*)', '', m.find('m:name', ns).text)})

    peptides = mascot_root.findall('m:hits/m:hit/m:protein/m:peptide', ns)

    rows_list = list()
    for pep in peptides:
        # Extract only first rank predictions.
        if pep.get('rank') == '1':
            row = dict()

            # Add scan title.
            row.update(
                {'title': pep.find('m:pep_scan_title', ns).text.strip()}
            )

            # Modify sequence with modifications, if any.
            seq = pep.find('m:pep_seq', ns).text
            mod_pos = pep.find('m:pep_var_mod_pos', ns).text
            if mod_pos is not None:
                mod_pos = re.search(
                    '\d\.(\d+)\.\d',
                    mod_pos
                ).group(1)
                seq = ''.join([aa if m == '0' else aa + '(' + mods[m] + ')'
                               for aa, m in zip(seq, mod_pos)])
            row.update({'sequence': seq})

            rows_list.append(row)

    return pd.DataFrame(rows_list)


def format_mgf_deepnovo(mgf_input: str, mgf_output: str):
    """Format MGF file for use with DeepNovoV2.

    Necessary spectrum parameters will be reordered to comply with DeepNovoV2
    convention. Other parameters will be discarded. Empty spectra will be
    discarded.

    Args:
        mgf_input (str): path to the input MGF file.
        mgf_output (str): path to the output MGF file.

    """
    key_order = ['title', 'pepmass', 'charge', 'scans', 'rtinseconds']
    with mgf.read(mgf_input, read_charges=False) as reader:
        for spectrum in reader:
            # Check if spectrum isn't emtpy.
            if spectrum['m/z array'].size:
                # Remove unnecessary parameters.
                to_remove = [c for c in spectrum['params'].keys() if c not in
                             key_order]
                for col in to_remove:
                    spectrum['params'].pop(col)
                # Append current spectrum to MGF output with correct params
                # order.
                mgf.write((spectrum,), mgf_output, key_order=key_order)
