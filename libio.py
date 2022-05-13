import argparse
import os
from pathlib import Path

def parse_arguments():
    """
    parse and check the command-line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--parameters", help="dat input parameter file")
    parser.add_argument("-v","--verbose", help="increase output verbosity",action="store_true")
    args = parser.parse_args()
    print(args.__dict__)
    # checks
    if not args.parameters:
        raise Exception("Parameter file is mandatory")
    if args.verbose:
        print("verbosity turned on")
    return args

def check_output_file(cleaned_pars):
    """If the output_filename already exists, block the execution"""
    output_path = Path(cleaned_pars["output_filename"])
    if output_path.is_file():
        raise Exception(f"Output file {cleaned_pars['output_filename']} already existing. Aborting")
    return

def system_parameters_setup(parfile):
    """Sets up the parameters."""
    pars = read_parfile(parfile)
    check_mandatory_parameters(pars)
    cleaned_pars = check_optional_parameters(pars)
    print(f"Cleaned Parameters {cleaned_pars}")
    check_output_file(cleaned_pars)
    return cleaned_pars

def read_parfile(parfile_string):
    # receives the name of the parameters file in input
    # parses it and gives back a dictionary
    ###################################################
    parfile = Path(parfile_string)
    print("reading parameters file ", parfile)
    if parfile.is_file == False:
        raise Exception("Parameter file not existing. Aborting.")
    parameters = {}
    with open(parfile, "r") as pars:
        for ln in pars:
            if ln.startswith("#") == False:
                split_list = ln.split()
                if len(split_list) != 2:
                    raise Exception(f"badly formatted parameter line\n{ln}")
                else:
                    par_name = split_list[0]
                    par_value = split_list[1]
                    parameters[par_name] = par_value
                    print("parameter ", par_name, " = ", par_value)
    return parameters


def check_mandatory_parameters(parameters):
    """Check the existence of mandatory parameters."""
    mandatory_keys = [
        "input_filename",
        "output_filename"
    ]

    observed_pars = parameters.keys()

    for par in mandatory_keys:
        if par not in observed_pars:
            raise Exception(f"missing parameter {par}")
    
    return


def check_optional_parameters(parameters):
    """Check the existence of optional parameters. Add their default value if absent"""
    optional_keys = {
        "max_binom" : ["integer", 100000],
    }
    
    observed_pars = parameters.keys()

    for optk in optional_keys.keys():
        if optk not in observed_pars:
            parameters[optk] = optional_keys[optk][1]
        else:
            if optional_keys[optk][0] == "integer":
                parameters[optk] = int(parameters[optk])
            if optional_keys[optk][0] == "float":
                parameters[optk] = float(parameters[optk])
        
    return parameters


def output_mappings(mapping_dict, mapping_order, output_filename):
    """
    Functions that outputs the mappings to a file
    
    Parameters
    ----------
    mapping_dict : dict
        dictionary of mappings
    
    mapping_order : list
        list of ordered mappings

    output_filename : str
    """
    
    header = "N\tmapping\ttrans_mapping\ths\thk\tsmap\tsmap_inf" + os.linesep
    with open(output_filename, "w") as wfile:
        # write header
        wfile.write(header)
        for ord_map in mapping_order:
            output_str = [] 
            for elem in mapping_dict[ord_map]:
                if isinstance(elem, float):
                    output_str.append(f"{elem:.6f}")
                else:
                    output_str.append(f"{elem}")
            wfile.write("\t".join(output_str) + os.linesep)
