import sys

from mfsolution import *
from analyse import *

class structtype():
    pass

if __name__ == "__main__":

    if (len(sys.argv) != 3) & (len(sys.argv) != 2 ) :
        raise ValueError('Invalid number of arguments. Usage: {} config_file.yaml output_path.yaml'.format(sys.argv[0]))

    print('Running solution code')

    config_file = sys.argv[1]

    mfsolution(reconfigure, config_file)

    if len(sys.argv) == 3:
        print('Running analysis code')

        output_path = sys.argv[2]
        analyse(fn_analyse, config_file,output_path)