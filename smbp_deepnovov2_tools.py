# -*- coding: utf-8 -*-
"""Create DeepNovoV2 input files from MGF spectra and Mascot annotation files.

This module contains tools to annotate MGF spectra using Mascot XML results and
create the features.csv file needed by DeepNovoV2.
"""

import os
import pandas as pd
import pyteomics.mgf as mgf
import re
from typing import List, Iterable
import xml.etree.ElementTree as ET


def extract_features(mgf_file: str):
    """Extract features from an MGF file to create DeepNovoV2 features.csv.

    Args:
        mgf_file (str): path to the MGF file.

    Returns:
        pandas.DataFrame: the features DataFrame.

    """
    with mgf.read(mgf_file, read_charges=False) as reader:
        rows_list = list()
        for spectrum in reader:
            params = spectrum['params']
            row = dict()
            # spec_group_id and scans are identical since we're not merging
            # spectra.
            row.update({'spec_group_id':  params['scans'],
                        'scans': params['scans']})

            # In Mascot MGF files, the PEPMASS line contains the precursor m/z,
            # intensity and charge: only m/z is selected.
            # Pyteomics returns charges as a list of integers with subtype
            # pyteomics.auxiliary.structures.Charge. The first charge is
            # selected and converted to integer.
            row.update({'m/z': params['pepmass'][0],
                        'z': int(params['charge'][0])})

            row.update({'rt_mean': params['rtinseconds']})
            row.update({'seq': 'NA'})

            # profile and feature area are set to default values set by the
            # authors.
            row.update({'profile': '0.0:1.0', 'feature area': '1.0'})

            rows_list.append(row)

    # Return DataFrame of features.
    return pd.DataFrame(rows_list)


def annotate_mgf(mgf_file: str, mascot_xml: str):
    """Annotate MGF file using Mascot XML results.

    annotate_mgf will annotate the MGF file using peptide sequences found in
    the Mascot XML results and write the resulting MGF file, appending
    '_annotated' to the original filename.

    Args:
        mgf_file (str): path to the MGF file.
        mascot_xml (str): path to the Mascot XML results.

    """
    def annotated_mgf_iter(mgf_reader: Iterable, mascot_xml: str):
        """Annotated MGF generator."""
        mascot_seq = extract_mascot_sequences(mascot_xml)
        for spectrum in mgf_reader:
            sequences = mascot_seq.loc[
                mascot_seq.title == spectrum['params']['title'], 'sequence'
            ].values
            for seq in sequences:
                spectrum['params']['seq'] = seq
                yield spectrum

    # Output MGF filename.
    output_mgf = '{0}_{2}{1}'.format(*os.path.splitext(mgf_file) +
                                     ('annotated',))
    with mgf.read(mgf_file, read_charges=False) as reader:
        mgf.write(annotated_mgf_iter(reader, mascot_xml), output_mgf)


def merge_mgf(input_file_list: List[str], output_file: str):
    """Merge a list of MGF files, re-numbering the scan IDs.

    Args:
        input_file_list (list): a list of paths to MGF files.
        output_file (str): path to the merged MGF file.

    """
    with open(output_file, mode="w") as output_handle:
        for index, input_file in enumerate(input_file_list):
            print("input_file = ", os.path.join(input_file))
            with open(input_file, mode="r") as input_handle:
                line = input_handle.readline()
                while line:
                    # If the current line is a scan ID, renumber it.
                    if "SCANS=" in line:
                        scan = re.split('=|\n', line)[1]
                        output_handle.write("SCANS=F{0}:{1}\n"
                                            .format(index, scan))
                    else:
                        output_handle.write(line)
                    line = input_handle.readline()


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
    convention. Other parameters will be discarded.

    Args:
        mgf_input (str): path to the input MGF file.
        mgf_output (str): path to the output MGF file.

    """
    key_order = ['title', 'pepmass', 'charge', 'scans', 'rtinseconds']
    with mgf.read(mgf_input, read_charges=False) as reader:
        for spectrum in reader:
            # Remove unnecessary parameters.
            to_remove = [c for c in spectrum['params'].keys() if c not in
                         key_order]
            for c in to_remove:
                spectrum['params'].pop(c)
            # Append current spectrum to MGF output with correct params order.
            mgf.write((spectrum,), mgf_output, key_order=key_order)
