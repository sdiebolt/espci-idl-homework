# -*- coding: utf-8 -*-
"""Merge SMBP MGF files and generate input files for DeepNovoV2."""

import glob
import os
import smbp_deepnovov2_tools as sdt


# Extract features.
mgf_files = glob.glob('smbp_data/*.dat.mgf')
for mgf_file in mgf_files:
    # Extract features for current MGF file.
    features = sdt.extract_features(
        mgf_file, os.path.splitext(mgf_file)[0]+'.xml'
    )
    # Save features to CSV, appending _features to the filename.
    features.to_csv(
        "{0}_{1}.csv".format(os.path.splitext(mgf_file)[0], 'features'),
        index=False
    )

# Merge MGF and features files.
features_files = glob.glob('smbp_data/*_features.csv')
sdt.merge_mgf(mgf_files, 'spectrum_smbp.mgf')
sdt.merge_features(features_files, 'features_smbp.mgf')

# Reformat MGF file to comply with DeepNovoV2.
sdt.format_mgf_deepnovo('spectrum_smbp.mgf', 'spectrum_smbp_formatted.mgf')
