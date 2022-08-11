"""Test the libentropy library."""

import os
from pathlib import Path

import numpy as np
import pytest

# import utils modules
from src.libs.libio import system_parameters_setup

from . import reference_data


@pytest.fixture
def example_parfile_measure():
    """Parameter file."""
    return Path(reference_data, "parameters_test_measure.dat")

@pytest.fixture
def example_parfile_optimize():
    """Parameter file for task optimize."""
    return Path(reference_data, "parameters_test_optimize.dat")

@pytest.fixture
def example_missing_parfile():
    """Non existing parameter file."""
    return Path(reference_data, "parameterf_test.dat")


@pytest.fixture
def example_incomplete_parfile():
    """Parameter file with no input_filename."""
    return Path(reference_data, "parameters_test_missing.dat")


@pytest.fixture
def example_existing_output_parfile():
    """Parameter file with already existing output_filename."""
    return Path(reference_data, "parameters_test_existing_output.dat")


def test_parameter_file_measure(example_parfile_measure):
    """Test correct functioning of system_parameter_setup."""
    expected_output_dict = {
        "input_filename": "input.csv",
        "output_filename": "output.csv",
        "max_binom": 2
    }

    observed_pars_dict = system_parameters_setup(example_parfile_measure,"measure")

    assert observed_pars_dict == expected_output_dict

def test_parameter_file_optimize(example_parfile_optimize):
    expected_output_dict = {
        "input_filename": "input.csv",
        "output_filename": "output.csv",
        "nsteps": 42,
        "ncg": 1
    }

    observed_pars_dict = system_parameters_setup(example_parfile_optimize,"optimize")

    assert observed_pars_dict == expected_output_dict
    

def test_missing_parameter_file(example_missing_parfile):
    """Check error if parameter file is missing."""
    with pytest.raises(Exception):
        system_parameters_setup(example_missing_parfile,"measure")


def test_incomplete_parameter_file(example_incomplete_parfile):
    """Check error if parameter file is missing."""
    with pytest.raises(Exception):
        system_parameters_setup(example_incomplete_parfile,"measure")


def test_existing_output(example_existing_output_parfile):
    """Check error if output filename already exists."""
    with pytest.raises(Exception):
        system_parameters_setup(example_existing_output_parfile, "measure")
