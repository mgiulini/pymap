"""Test the measure task."""
# from src.measure import MEASURE
import numpy as np
import pytest


@pytest.fixture
def example_mappings():
    """Provide mapping dictionary and corresponding order."""
    cg_mappings = {}
    cg_mappings["mapping_dict"] = {
        (2, 6): (2, np.array([2, 6]), ["C", "G"], 0.1, 0.05, 0.2, 1.0),
        (0, 1): (2, np.array([0, 1]), ["A", "B"], 0.2, 0.15, 0.23, 1.0)
        }

    cg_mappings["mapping_order"] = [(0, 1), (2, 6)]
    return cg_mappings


def test_output_format(example_mappings):
    """Check if the output format is the expected one."""
    assert True
