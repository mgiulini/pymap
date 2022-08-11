"""Library to perform input-output tasks."""
import argparse
from pathlib import Path

TASKS = ["measure", "optimize"]

def parse_arguments():
    """Parse and check the command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--task", help="task to be executed")
    parser.add_argument("-p", "--parameters", help="dat input parameter file")
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true"
    )
    args = parser.parse_args()
    print(args.__dict__)
    # checks
    if not args.task:
        raise Exception("Task is mandatory")
    if args.task not in TASKS:
        raise Exception(f"Invalid task, should be in {' '.join(TASKS)}")
    if not args.parameters:
        raise Exception("Parameter file is mandatory")
    if args.verbose:
        print("verbosity turned on")
    return args


def check_output_file(cleaned_pars):
    """If the output_filename already exists, block the execution."""
    output_fl = cleaned_pars["output_filename"]
    output_path = Path(output_fl)
    print(f"checking output path {output_path}")
    if output_path.is_file():
        raise Exception(f"Output file {output_fl} already existing. Aborting")
    return


def system_parameters_setup(parfile, task):
    """Set up the parameters."""
    pars = read_parfile(parfile)
    check_mandatory_parameters(pars)
    cleaned_pars = check_optional_parameters(pars, task)
    print(f"Cleaned Parameters {cleaned_pars}")
    check_output_file(cleaned_pars)
    return cleaned_pars


def read_parfile(parfile_string):
    """Read parameter file."""
    parfile = Path(parfile_string)
    print("reading parameters file ", parfile)
    if parfile.is_file is False:
        raise Exception("Parameter file not existing. Aborting.")
    parameters = {}
    with open(parfile, "r") as pars:
        for ln in pars:
            if ln.startswith("#") is False:
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


def check_optional_parameters(parameters, task):
    """
    Check the existence of task-specific optional parameters.
    Set to default if absent.
    """
    optional_keys = {
        "measure" : {"max_binom": ["integer", 100000]},
        "optimize" : {
            "nsteps" : ["integer", 100],
            "ncg" : ["integer", 1] # default (not so useful) choice is 1
        }
    }

    observed_pars = parameters.keys()

    for optk in optional_keys[task].keys():
        if optk not in observed_pars:
            # paramater not present. set to default
            parameters[optk] = optional_keys[task][optk][1]
        else:
            # parameter is present, convert to the desired type
            if optional_keys[task][optk][0] == "integer":
                parameters[optk] = int(parameters[optk])
            if optional_keys[task][optk][0] == "float":
                parameters[optk] = float(parameters[optk])
    return parameters
