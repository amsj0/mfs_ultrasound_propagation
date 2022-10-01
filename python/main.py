import sys

from mfsolution import *
from analyse import *
from tseries import *

if __name__ == "__main__":

    #Store ST
    
    if (len(sys.argv) != 4) & (len(sys.argv) != 3) & (len(sys.argv) != 2 ) :
        raise ValueError('Invalid number of arguments. Usage: {} inputfile.yaml output_path.yaml config.yaml'.format(sys.argv[0]))

    print('Running solution code')

    input_file = sys.argv[1]

    ST = mfsolution(__name__,input_file)

    if len(sys.argv) > 2:
        print('Running analysis code')

        output_path = sys.argv[2]
        analyse(ST,input_file,output_path)

    if len(sys.argv) == 4:
        print('Running tseries code')

        config_file = sys.argv[3]
        tseries(config_file, output_path)