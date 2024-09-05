import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import upsetplot

#######################
### INPUT FUNCTIONS ###
#######################

def parse_args():
    """ Parse commandline arguments """
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--in_config", default=None, type=str,
                        required=True, help="The input csv")

    return parser.parse_args()

def parse_input_config(in_config):
    """ This will convert the input config into a dictionary """

    in_config_dict = {}

    with open(in_config, 'r') as f_in:
        for line in f_in.readlines():
            prefix, file_path = line.strip().split(',')
            in_config_dict[prefix] = file_path

    return in_config_dict

##########################
### ANALYSIS FUNCTIONS ###
##########################

def analyze_data(bc_file_dict):
    """ This will process the data passed in and generate a variety of plots
        and figures to be used to assess any number of barcode lists in order
        to determine how similar the lists are to each other

    Args:
        whitelist (list): The barcode whitelist
        bc_file_dict (dict): The barcode files to analyse with associated prefixes

    Returns:
        None
    """

    dict_for_upset = {}
    for bc_file_tag, bc_file in bc_file_dict.items():
        prefix = bc_file_tag
        bc_df = pd.read_csv(bc_file, header=None, names=['bc'])
        dict_for_upset[prefix] = bc_df.bc.tolist()

    # Upset Plot
    make_upset_plot(dict_for_upset)

def make_upset_plot(counts_dict):
    """ Makes an upset plot for the data

    Args:
        counts_dict (dict): The input for the upset plot

    Returns:
        None
    """
    upsetplot.plot(upsetplot.from_contents(counts_dict), show_counts = True)
    plt.suptitle("5' Barcode Comparison")
    plt.savefig('5_prime_bcs.pdf')
    plt.close()

#######################
### MAIN SUBROUTINE ###
#######################

def main():
    args = parse_args()

    bc_file_dict = parse_input_config(args.in_config) 

    analyze_data(bc_file_dict)

if __name__ == '__main__':
    main()
