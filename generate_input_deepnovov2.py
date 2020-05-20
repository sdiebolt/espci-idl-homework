# -*- coding: utf-8 -*-
"""Merge SMBP MGF files and generate input files for DeepNovoV2."""

import glob
import os
import smbp_deepnovov2_tools as sdt


# Extract features.
mgf_files = glob.glob('smbp_data/*.dat.mgf')
for mgf_file in mgf_files:
    print('Extracting features from: {}'.format(mgf_file))

    mgf_file_formatted = '{0}_formatted{1}'.format(*os.path.splitext(mgf_file))
    mascot_file = os.path.splitext(mgf_file)[0]+'.xml'
    features_file = (
        "{0}_features.csv"
        .format(os.path.splitext(mgf_file_formatted)[0])
    )

    # Reformat MGF file to comply with DeepNovoV2.
    sdt.format_mgf_deepnovo(mgf_file, mgf_file_formatted)
    print('Formatted MGF: {}'.format(mgf_file_formatted))

    # Extract features for current MGF file.
    features = sdt.extract_features(
        mgf_file_formatted, mascot_file
    )

    # Save features to CSV, appending _features to the filename.
    features.to_csv(features_file, index=False)
    print('Extracted features: {}\n'.format(features_file))

# Merge MGF and features files.
mgf_formatted_files = sorted(glob.glob('smbp_data/*_formatted.mgf'))
features_files = sorted(glob.glob('smbp_data/*_features.csv'))
sdt.merge_mgf(mgf_formatted_files, 'spectrum_smbp.mgf')
sdt.merge_features(features_files, 'features_smbp.csv')

sdt.partition_feature_file_nodup('features_smbp.csv', prob=[0.8, 0.1, 0.1])
