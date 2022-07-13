import sys

from mfsolution import *
from analyse import *
from tseries import *

class structtype():
    pass

if __name__ == "__main__":

    if (len(sys.argv) != 3) & (len(sys.argv) != 3) & (len(sys.argv) != 2 ) :
        raise ValueError('Invalid number of arguments. Usage: {} inputfile.yaml output_path.yaml config.yaml'.format(sys.argv[0]))

    print('Running solution code')

    input_file = sys.argv[1]

    mfsolution(reconfigure, input_file)

    if len(sys.argv) > 2:
        print('Running analysis code')

        output_path = sys.argv[2]
        analyse(fn_analyse, input_file,output_path)

    if len(sys.argv) == 4:
        print('Running tseries code')

        config_file = sys.argv[3]
        tseries(pre_config, set_domain_plot, create_matrix, plot_save_table, config_file, output_path)