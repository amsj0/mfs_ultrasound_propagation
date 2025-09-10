import sys

from mfsolution import *
from analyse_concurrent import *
from run_tseries_parallel import *

if __name__ == "__main__":

    #Store ST
    
    if (len(sys.argv) != 4) & (len(sys.argv) != 3) & (len(sys.argv) != 2 ) :
        raise ValueError('Invalid number of arguments. Usage: {} inputfile.yaml output_path.yaml config.yaml'.format(sys.argv[0]))

    print('Running solution code')

    input_file = sys.argv[1]
    output_path = ''

    ST = mfsolution(input_file,output_path,__name__)

    if len(sys.argv) > 2:
        print('Running analysis code')

        output_path = sys.argv[2]
        ST = analyse_concurrent(
            store=ST,
            input_file=input_file,
            output_path=output_path,
            task_name=__name__
        )

    if len(sys.argv) == 4:
        print('Running tseries code')

        config_file = sys.argv[3]
        run_tseries_parallel(
            store=ST,
            config_file=config_file,
            output_path=output_path,
            task_name=__name__
        )