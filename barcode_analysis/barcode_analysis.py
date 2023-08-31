import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import upsetplot

MIN_Q_SCORE = 15
EXP_COUNT = 1000

#######################
### INPUT FUNCTIONS ###
#######################

def parse_args():
    """ Parse commandline arguments """
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--in_config", default=None, type=str,
                        required=True, help="The input csv")

    parser.add_argument("-w", "--whitelist", default=None, type=str,
                        required=True, help="Barcode whitelist")

    return parser.parse_args()

def parse_input_config(in_config):
    """ This will convert the input config into a dictionary """

    in_config_dict = {}

    with open(in_config, 'r') as f_in:
        for line in f_in.readlines():
            prefix, file_path = line.strip().split(',')
            in_config_dict[prefix] = file_path

    return in_config_dict

def read_whitelist(whitelist_file):
    """ Reads in the whitelist file into a list

    Args:
        whitelist_file (str): The path to the barcode whitelist file

    Returns:
        whitelist_list (list): The barcode whitelist
    """
    whitelist_list = []

    with open(whitelist_file) as f:
        for line in f:
            whitelist_list.append(line.strip())

    return set(whitelist_list)


##########################
### ANALYSIS FUNCTIONS ###
##########################

def analyze_data(whitelist, bc_file_dict):
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
        bc_df = pd.read_csv(bc_file)
        
        confident_bcs = get_confident_bcs(whitelist, MIN_Q_SCORE, bc_df)

        # Barcode Rank Plot
        threshold = bc_quantile(confident_bcs.bc_count, EXP_COUNT)
        make_knee_plot(prefix, threshold, confident_bcs.bc_count)

        # Make inputs for Upset Plot
        dict_for_upset[prefix] = make_upset_plot_input(confident_bcs, threshold)

    # Upset Plot
    make_upset_plot(dict_for_upset)

def bc_quantile(counts, exp_count):
    count = np.sort(counts)[::-1][:exp_count]
    return np.quantile(count, 0.95) / 20

def get_confident_bcs(whitelist, min_q_score, bc_df):
    """ Filter the estimated barcodes using the min_q_score and the whitelist

    Args:
        whitelist (list): The list of whitelisted barcodes
        min_q_score (int): The minimum q score for a barcode to be considered 'correct'
        bc_df (dataframe): 

    Returns:
        filtered_bc_df (dataframe)
    """

    filtered_bcs = bc_df[bc_df.bc_min_q >= min_q_score]
    filtered_bcs = pd.DataFrame(filtered_bcs.bc.value_counts().reset_index())
    filtered_bcs.columns = ['bc', 'bc_count']
    filtered_bcs = filtered_bcs[filtered_bcs.bc.isin(whitelist)]

    return filtered_bcs

def make_knee_plot(prefix, threshold, counts):
    """ Makes a knee plot of the data

    Args:
        prefix (str): The prefix for the output file name
        threshold (int): The counts threshold used for determining high
            confidence barcodes
        counts (pd.Series): The barcode counts

    Returns:
        None
    """
    fig, ax = plt.subplots(1, 1, figsize = (8, 5))

    ax.loglog(list(counts), marker = "o", linestyle = "")

    ax.set_xlabel('Barcodes')
    ax.set_ylabel('High Confident Counts')
    ax.set_title('Knee Plot')
    ax.axhline(y=threshold, color='r', linestyle='--', label = 'Cell Calling Threshold')
    ax.legend()

    fig.savefig(prefix + '_knee_plot.png')
    plt.close()

def make_upset_plot_input(filtered_bc, threshold):
    """ Converts the dataframe so its able to be used by UpSetPlot

    Args:
        filtered_bc (dataframe): The dataframe containing filtered bcs
        threshold (int): The counts threshold used for determining high
            confidence counts

    Returns:
        threshold_bc (list): The barcodes that are over the threshold
    """
    threshold_counts = filtered_bc.bc_count >= threshold
    threshold_bc = filtered_bc[threshold_counts].bc.tolist()

    return threshold_bc

def make_upset_plot(counts_dict):
    """ Makes an upset plot for the data

    Args:
        counts_dict (dict): The input for the upset plot

    Returns:
        None
    """
    upsetplot.plot(upsetplot.from_contents(counts_dict), show_counts = True)
    plt.savefig('upset.png')
    plt.close()

#######################
### MAIN SUBROUTINE ###
#######################

def main():
    args = parse_args()

    bc_file_dict = parse_input_config(args.in_config) 
    whitelist = read_whitelist(args.whitelist)

    analyze_data(whitelist,
                 bc_file_dict)

if __name__ == '__main__':
    main()
