# -*- coding: utf-8 -*-
"""Merge SMBP MGF files and generate input files for DeepNovoV2."""

import glob
import os
import smbp_deepnovov2_tools as sdt


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
