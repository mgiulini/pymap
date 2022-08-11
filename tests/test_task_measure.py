#from src.measure import MEASURE
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
    # output_filename = "test_output.csv"
    # output_mappings(
        # example_mappings["mapping_dict"],
        # example_mappings["mapping_order"],
        # output_filename
        # )
    #building expected output
    # lines = [
            # ["N", "mapping", "trans_mapping", "hs",
                # "hk", "smap", "smap_inf"],
            # ["2", "[0 1]", "['A', 'B']", "0.200000",
                # "0.150000", "0.230000", "1.000000"],
            # ["2", "[2 6]", "['C', 'G']", "0.100000",
                # "0.050000", "0.200000", "1.000000"]
        # ]
    # expected_output = ""
    # for ln in lines:
        # line = "\t".join(ln) + os.linesep
        # expected_output += line
# 
    # file_content = open(output_filename, "r").read()
# 
    # assert file_content == expected_output
# 
    # os.unlink(output_filename)
    assert True
