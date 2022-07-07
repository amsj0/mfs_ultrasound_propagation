import sys

from mfsolution import *
from analyse import *

class structtype():
    pass

if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise ValueError('Invalid number of arguments. Usage: {} config_file.yaml'.format(sys.argv[0]))

    yaml_path = sys.argv[1]

    mfsolution(reconfigure, yaml_path)
    analyse(fn_analyse, yaml_path)