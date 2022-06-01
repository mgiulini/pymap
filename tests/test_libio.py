"""Test the libentropy library."""

import os
from pathlib import Path
import numpy as np
import pytest
# import utils modules
from pymap.libio import (
    output_mappings,
    system_parameters_setup,
    parse_arguments,
)

from . import reference_data

@pytest.fixture
def example_parfile():
    """example parameter file."""
    return Path(reference_data, "parameters_test.dat")

@pytest.fixture
def example_missing_parfile():
    """non existing parameter file."""
    return Path(reference_data, "parameterf_test.dat")

@pytest.fixture
def example_incomplete_parfile():
    """parameter file with no input_filename."""
    return Path(reference_data, "parameters_test_missing.dat")

@pytest.fixture
def example_existing_output_parfile():
    """parameter file with already existing output_filename."""
    return Path(reference_data, "parameters_test_existing_output.dat")

@pytest.fixture
def example_mappings():
    """Example mapping dictionary and corresponding order."""
    cg_mappings = {}
    cg_mappings["mapping_dict"] = {
        (2,6) : (2, np.array([2,6]), ["C", "G"], 0.1, 0.05, 0.2, 1.0),
        (0,1) : (2, np.array([0,1]), ["A", "B"], 0.2, 0.15, 0.23, 1.0)
        }
    
    cg_mappings["mapping_order"] = [(0,1),(2,6)]
    return cg_mappings

def test_parameter_file(example_parfile):
    """Test correct functioning of system_parameter_setup."""
    expected_output_dict = {
        "input_filename" : "input.csv",
        "output_filename" : "output.csv",
        "max_binom" : 2
    }

    observed_pars_dict = system_parameters_setup(example_parfile)

    assert observed_pars_dict == expected_output_dict

def test_missing_parameter_file(example_missing_parfile):
    """Check error if parameter file is missing."""
    with pytest.raises(Exception):
        system_parameters_setup(example_missing_parfile)

def test_incomplete_parameter_file(example_incomplete_parfile):
    """Check error if parameter file is missing."""
    with pytest.raises(Exception):
        system_parameters_setup(example_missing_parfile)

def test_existing_output(example_existing_output_parfile):
    """Check error if output filename already exists"""
    with pytest.raises(Exception):
        system_parameters_setup(example_existing_output_parfile)

def test_output_format(example_mappings):
    """Check if the output format is the expected one."""
    output_filename = "test_output.csv"
    output_mappings(
        example_mappings["mapping_dict"],
        example_mappings["mapping_order"],
        output_filename
        )
    # building expected output
    lines = [
        "N", "mapping", "trans_mapping", "hs", "hk", "smap", "smap_inf",
        "2", "[0 1]", "['A', 'B']", "0.200000", "0.150000", "0.230000", "1.000000",
        "2", "[2 6]", "['C', 'G']", "0.100000", "0.050000", "0.200000", "1.000000"
        ]
    expected_output = ""
    for ln in lines:
        line = "\t".join(ln) + os.linesep
        expected_output += line

    #expected_output = "N\tmapping\ttrans_mapping\ths\thk\tsmap\tsmap_inf" + os.linesep
    #expected_output += "2\t[0 1]\t['A', 'B']\t0.200000\t0.150000\t0.230000\t1.000000" + os.linesep
    #expected_output += "2\t[2 6]\t['C', 'G']\t0.100000\t0.050000\t0.200000\t1.000000" + os.linesep

    file_content = open(output_filename, "r").read()

    assert file_content == expected_output

    os.unlink(output_filename)