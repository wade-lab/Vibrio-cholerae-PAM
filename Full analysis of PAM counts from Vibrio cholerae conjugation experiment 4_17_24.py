import tkinter as tk
from tkinter import filedialog
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

def select_file(prompt):
    """
    Opens a file dialog to select a file based on a given prompt.
    
    Args:
        prompt (str): The prompt displayed in the file dialog.
    
    Returns:
        str: The path to the selected file.
    """
    root = tk.Tk()  # Create a new Tkinter root window.
    root.withdraw()  # Hide the root window immediately.
    file_path = filedialog.askopenfilename(
        title=prompt,
        filetypes=(("FASTQ files", "*.fastq"), ("All files", "*.*"))  # Set filter for FASTQ files and all files.
    )
    return file_path

def count_PAMs(filename, possible_PAMs, seq_start='AGGAAACAGCTATGACCATGATTACGCCAAGCTT'): #expected sequence at the start of each read
    """
    Counts occurrences of PAM sequences in a file that start after a specific sequence.
    
    Args:
        filename (str): Path to the file containing sequences.
        possible_PAMs (list): List of PAM sequences to count.
        seq_start (str): The sequence that precedes the PAM sequences.
        
    Returns:
        list: Counts of each PAM in the order they appear in possible_PAMs.
    """
    with open(filename, 'r') as file:
        lines = file.read().splitlines()

    # Extracting the PAM sequences that follow directly after seq_start in each line. Includes 2 nt before the actual PAM to account for the empty vector control sequence
    list_of_extended_PAMs = [line[:39][34:] for line in lines if line.startswith(seq_start)]
    extended_PAM_count = Counter(list_of_extended_PAMs)  # Count each unique PAM sequence.

    return [extended_PAM_count[pam] for pam in possible_PAMs]

def normalize_data(file_details, possible_PAMs):
    """
    Normalizes raw PAM counts based on comparison between 'active' and 'inactive' CRISPR conditions.

    Args:
        file_details (dict): A dictionary mapping descriptive labels to file paths.
        possible_PAMs (list): A list of all possible PAMs.

    Returns:
        tuple: A tuple containing dictionaries for normalized counts, standard deviations, and raw counts.
    """
    # Load counts for all files, described by 'file_details'.
    pam_counts = {description: count_PAMs(filename, possible_PAMs) for description, filename in file_details.items()}
    normalized_counts = {}

    # Normalize 'active' counts by their corresponding 'inactive' counts.
    for description in file_details.keys():
        if 'active' in description and 'inactive' not in description:
            inactive_description = description.replace('active', 'inactive')
            inactive_counts = np.array(pam_counts[inactive_description])
            active_counts = np.array(pam_counts[description])
            normalized_counts[description] = active_counts / inactive_counts

    # Normalize further by the 'GCATG' PAM value.
    for description, counts in normalized_counts.items():
        gcatg_index = possible_PAMs.index('GCATG')
        gcatg_value = counts[gcatg_index]
        normalized_counts[description] = counts / gcatg_value

    # Compute the mean and standard deviation for each PAM.
    average_normalized_counts = {}
    std_devs = {}
    non_gcatg_indices = list(range(len(possible_PAMs)))
    non_gcatg_indices.remove(gcatg_index)

    for pam_index in non_gcatg_indices:
        for description, counts in normalized_counts.items():
            base_description = "Spacer " + description.split(' ')[1].replace(',', '')
            if base_description not in average_normalized_counts:
                average_normalized_counts[base_description] = [[] for _ in non_gcatg_indices]
                std_devs[base_description] = [[] for _ in non_gcatg_indices]
            average_normalized_counts[base_description][pam_index].append(counts[pam_index])

    # Calculate averages and standard deviations.
    for key in average_normalized_counts:
        averages = []
        standard_devs = []
        for pam in average_normalized_counts[key]:
            averages.append(np.mean(pam))
            standard_devs.append(np.std(pam))
        average_normalized_counts[key] = averages
        std_devs[key] = standard_devs

    return average_normalized_counts, std_devs, pam_counts

def write_output(average_normalized_counts, possible_PAMs, output_file):
    """
    Writes the normalized conjugation efficiencies to an output file.

    Args:
        average_normalized_counts (dict): A dictionary of normalized PAM counts.
        possible_PAMs (list): A list of PAM sequences.
        output_file (str): Path to the output file.
    """
    non_gcatg_pams = [pam[2:5] for pam in possible_PAMs if pam != "GCATG"]
    with open(output_file, 'w') as f:
        f.write('PAM\tSpacer 4\tSpacer 21\n')
        for i, pam in enumerate(non_gcatg_pams):
            f.write(f'{pam}\t{average_normalized_counts["Spacer 4"][i]}\t{average_normalized_counts["Spacer 21"][i]}\n')

def write_output_raw_counts(raw_PAM_counts, possible_PAMs, output_file_raw_counts):
    """
    Writes raw PAM counts to a designated output file.

    Args:
        raw_PAM_counts (dict): Raw counts of each PAM sequence.
        possible_PAMs (list): A list of all possible PAMs.
        output_file_raw_counts (str): Path to the raw counts output file.
    """
    with open(output_file_raw_counts, 'w') as f:
        f.write('\t')
        for PAM in possible_PAMs:
            if PAM[0:2] == 'AT':
                f.write(PAM[2:5] + '\t')
            elif PAM == 'GCATG':
                f.write('Empty vector control\t')
        f.write('\n')
        for label in raw_PAM_counts:
            f.write(label + '\t')
            for read_coverage in raw_PAM_counts[label]:
                f.write(str(read_coverage) + '\t')
            f.write('\n')

def create_scatter_plot(data, std_devs, output_path, pam_sequences):
    """
    Creates and saves a scatter plot of the normalized data with error bars.

    Args:
        data (dict): Normalized data for plotting.
        std_devs (dict): Standard deviations for error bars.
        output_path (str): Path to save the plot image.
        pam_sequences (list): List of PAM sequences for color coding.
    """
    plt.figure()
    x_data = data['Spacer 4']
    y_data = data['Spacer 21']
    x_errors = std_devs['Spacer 4']
    y_errors = std_devs['Spacer 21']

    # Define PAM groups and their respective colors.
    pam_groups = {
        'hot1': ['ATAAC', 'ATAAT'],
        'hot2': ['ATATT', 'ATATC'],
        'medium': ['ATTAC', 'ATCAT', 'ATCAC', 'ATGAT', 'ATTAT', 'ATGAC', 'ATAGC', 'ATAAG', 'ATAGT'],
        'cool1': ['ATAAA'],
        'cool2': ['ATAGG'],
    }
    colors = {
        'hot1': 'red',
        'hot2': 'orange',
        'medium': 'olive',
        'cool1': 'darkcyan',
        'cool2': 'slateblue',
        'other': 'gray'  # Default color for ungrouped PAMs
    }

    # Assign colors based on PAM group membership.
    point_colors = []
    for pam in pam_sequences:
        found = False
        for group, members in pam_groups.items():
            if pam in members:
                point_colors.append(colors[group])
                found = True
                break
        if not found:
            point_colors.append(colors['other'])
    
    # Plot each point with corresponding error bars and colors.
    for i in range(len(x_data)):
        plt.errorbar(x_data[i], y_data[i], xerr=x_errors[i], yerr=y_errors[i], fmt='o', color=point_colors[i], ecolor=point_colors[i], elinewidth=1, capsize=2)

    plt.xlabel('Normalized Conjugation Efficiency\nSpacer 4')
    plt.ylabel('Normalized Conjugation Efficiency\nSpacer 21')
    plt.gca().set_aspect('equal', adjustable='datalim') # Set aspect of the plot to be equal

    tick_interval = 0.25  # Adjust as needed for your data range and scale
    plt.gca().xaxis.set_major_locator(MultipleLocator(tick_interval))
    plt.gca().yaxis.set_major_locator(MultipleLocator(tick_interval))

    plt.subplots_adjust(bottom=0.15)  # Adjust layout to prevent label cutoff
    plt.savefig(output_path.replace('.txt', '_scatter.png'))
    plt.show()


def main():
    """
    Main function to orchestrate file selection, data processing, and output generation.
    """
    DNA_bases = ['A', 'C', 'G', 'T']
    possible_PAMs = ['AT' + a + b + c for a in DNA_bases for b in DNA_bases for c in DNA_bases] + ['GCATG']
    PAM_sequences = ['AT' + a + b + c for a in DNA_bases for b in DNA_bases for c in DNA_bases]
    prompts = [
        "Spacer 4, CRISPR-inactive, replicate 1",
        "Spacer 4, CRISPR-active, replicate 1",
        "Spacer 21, CRISPR-inactive, replicate 1",
        "Spacer 21, CRISPR-active, replicate 1",
        "Spacer 4, CRISPR-inactive, replicate 2",
        "Spacer 4, CRISPR-active, replicate 2",
        "Spacer 21, CRISPR-inactive, replicate 2",
        "Spacer 21, CRISPR-active, replicate 2"
    ]

    file_details = {prompt: select_file(prompt) for prompt in prompts}
    output_file_raw_counts = filedialog.asksaveasfilename(title='Save output file for raw PAM counts', filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    output_file = filedialog.asksaveasfilename(title='Save output file for normalized conjugation efficiencies', filetypes=[("Text files", "*.txt"), ("All files", "*.*")])

    if not output_file.endswith('.txt'):
        output_file += '.txt'
    if not output_file_raw_counts.endswith('.txt'):
        output_file_raw_counts += '.txt'

    average_normalized_counts, std_devs, raw_PAM_counts = normalize_data(file_details, possible_PAMs)
    write_output(average_normalized_counts, possible_PAMs, output_file)
    write_output_raw_counts(raw_PAM_counts, possible_PAMs, output_file_raw_counts)
    create_scatter_plot(average_normalized_counts, std_devs, output_file, PAM_sequences)

if __name__ == "__main__":
    main()
