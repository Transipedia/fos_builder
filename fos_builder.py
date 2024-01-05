#!/usr/bin/env python3

"""
Aims: build fos.txt file.It contains ordered list of samples proceeded by Reindeer, with or
whithout kmers found information picked on outpur of bcalm.
kmers found will be used to normalize Reindeer query output.
Input:
  - fof_unitigs.txt generated to build Reindeer index (required).
  - bcalm log file (optional).
    they must contain 'kmers found' value (expected bcalm log file format: <sample>_bcalm.log).
output:
    stdin or file whith lines as format:
    <sample name>   <number of kmers found by bcalm>
"""


import os
import sys
import argparse
import json
import re


__appname__   = "bcalm-tools"
__shortdesc__ = "Short description."
__licence__   = "none"
__version__   = "0.4.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


SAMPLE_NAME = ('Sample')
### targeted fields of multiqc files
LENGTH      = ('avg_sequence_length', 'FastQC_mqc-generalstats-fastqc-avg_sequence_length')
READ_COUNT  = ('FastQC_mqc-generalstats-fastqc-total_sequences', 'Total Sequences')
FILE        = 'multiqc_fastqc.txt'      # or 'multiqc_data.json'


def __main():
    """ Function doc """
    args = __usage()
    ### find provided multiqc find
    is_found, multiqc_file = get_filepath(args.multiqc_dir, FILE)
    if not is_found:
        sys.exit(f"{Col.RED}Error: {multiqc_file}")
    ### get matched samples and associated extra info
    samples, extra_info = get_samples(args, args.fof, multiqc_file)
    ### Output
    output_results(samples, extra_info, args.output)


def get_samples(args, fof, multiqc_file=None, bcalm_log_dir=None, prefix='', suffix='_bcalm.log'):
    """ Function doc """
    extra_info = None       # only if multiqc if provided
    ### get unitigs in right order
    try:
        with open(fof) as file:
            unitigs_path = file.readlines()
    except FileNotFoundError:
        print(f"Error: file {fof} not found.")
        sys.exit(1)
    ### extract sample names
    samples = []
    for sample in unitigs_path:
        samples.append(os.path.basename(sample)[:-12].rstrip())
    # ~ ### if bcalm_log argument is provided, add number of kmers
    # ~ if bcalm_log_dir:
        # ~ samples = __get_kmers_found(bcalm_log_dir, samples, prefix, suffix)
    ### if a multiqc file is provided, add reads information to be abble to normalize results
    if multiqc_file:
        samples, extra_info = __get_multiqc_info(args, multiqc_file, samples)
        return samples, extra_info
    return samples, extra_info


def __get_multiqc_info(args, multiqc_file, samples):
    """ get number and lenght reads found in samples of study (single or paired-end) """
    if not os.path.isfile(multiqc_file):
        print(f"Error: file {multiqc_file} not found. Read infos will not be reported.", file=sys.stderr)
        return samples
    with open(multiqc_file) as file:
        type_mqc = None
        if os.path.splitext(multiqc_file)[1].lower() == '.json':
            data = json.load(file)
            samples, extra_info = __json_info(args, samples, data)
        else:
            data = file.readlines()
            samples, extra_info = __text_info(args, samples, data)
    return samples, extra_info


def __text_info(args, samples, data):
    """ Function doc """
    total_read_num = 0
    total_kmer_num = 0
    header = data[0].strip().split('\t')
    ### find indexes of sample name, number of reads and read length.
    idx_length, idx_read_count = None, None
    idx_sample = header.index(SAMPLE_NAME)
    for length in LENGTH:
        if length in header:
            idx_length = header.index(length)
            break
    for read_count in READ_COUNT:
        if read_count in header:
            idx_read_count = header.index(read_count)
            break
    # ~ print('SAMPLE:', idx_sample,'LENGTH:', idx_length, 'COUNT:', idx_read_count, file=sys.stderr)
    ### aggregate infos
    samples_with_multiqc_infos = []
    for sample in samples:
        kmer_num = 0
        for line in data[1:]:
            ### When samples are in single mode
            infos = line.split('\t')
            if infos[idx_sample] == sample:
                read_len = int(float(infos[idx_length]))
                read_num = int(float(infos[idx_read_count]))
                kmer_num = (read_len - args.kmer_len + 1) * read_num
                total_read_num += read_num
            ### When samples are in paired-end mode
            if infos[idx_sample] == sample + '_1':
                read_len = int(float(infos[idx_length]))
                read_num = int(float(infos[idx_read_count]))
                kmer_num = (read_len - args.kmer_len + 1) * read_num * 2      # * 2 because of paired-end
                total_read_num += read_num
        ### if missing only one multiqc info, none of the mutliqc infos will be reported
        if kmer_num == 0:
            sys.exit(f"{Col.RED}Error: no multiqc info found for sample '{sample}'.")
        samples_with_multiqc_infos.append(f"{sample}\t{kmer_num}")
        total_kmer_num += kmer_num
    return samples_with_multiqc_infos, (total_read_num, total_kmer_num)


def __json_info(args, samples, data):
    """ Function doc """
    samples_with_multiqc_infos = []
    total_read_num = 0
    total_kmer_num = 0
    for sample in samples:
        kmer_num = 0
        for fastq_infos in data["report_general_stats_data"]:
            for fastq, infos in fastq_infos.items():
                ### When samples are in single mode
                if fastq == sample:
                    read_len = int(infos['avg_sequence_length'])
                    read_num = int(infos['total_sequences'])
                    kmer_num = (read_len - args.kmer_len + 1) * read_num
                    total_read_num += read_num
                ### When samples are in paired-end mode
                elif fastq == sample + '_1':
                    read_len = int(infos['avg_sequence_length'])
                    read_num = int(infos['total_sequences'])
                    kmer_num = (read_len - args.kmer_len + 1) * read_num * 2   # * 2 because of paired-end
                    total_read_num += read_num
        ### if missing only one multiqc info, none of the mutliqc infos will be reported
        if kmer_num == 0:
            sys.exit(f"{Col.RED}Error: no multiqc info found for sample {sample}.{Col.END}")
        samples_with_multiqc_infos.append(f"{sample}\t{kmer_num}")
        total_kmer_num += kmer_num
    return samples_with_multiqc_infos, (total_read_num, total_kmer_num)


def __get_kmers_found(bcalm_log_dir, samples, prefix, suffix):
    """ Function doc """

    bcalm_log_files = os.listdir(bcalm_log_dir)

    ### There must be at least the same count of bcalm log files as samples
    if len(bcalm_log_files) < len(samples):
        print("The number of bcalm log files must at least equal as the number of samples",
              f"({len(bcalm_log_files)} log files vs {len(samples)} samples).", file=sys.stderr)
        print("kmers found will not be reported.", file=sys.stderr)
        return samples

    ### associate kmers found count at sample
    samples_with_kmers_count = []
    for sample in samples:
        bcalm_log_file = os.path.join(bcalm_log_dir, prefix + sample + suffix)
        try:
            with open(bcalm_log_file) as file:
                for line in file:
                    if "kmers found" in line:
                        kmers_found = line.split()[6]
                        samples_with_kmers_count.append(f"{sample}\t{kmers_found}")
                        break
        except FileNotFoundError:
            print(f"{col.RED}Error: file {bcalm_log_file} not found.\nkmers found will not be reported.{Col.END}",
                  file=sys.stderr)
            return samples
    return samples_with_kmers_count


def output_results(samples, extra_info, output=None):
    """ Function doc """
    if output:
        with open(output, 'w') as file:
            file.write("{}\n".format('\n'.join(samples)))
    else:
        print(*samples, sep='\n')
    ### Print extra info
    if extra_info:
        print(f"\n{Col.GREEN}Extra info\n----------\n  Total reads: {extra_info[0]}\n  Total kmers: {extra_info[1]}.{Col.END}", file=sys.stderr)


def get_filepath(search_dir, file_name):
    """
    input:
        - search_dir: the search directory, it can be the file himself
        - file_name: the file to found
    output:
        - is_found: True or False
        - filepath: full path of the file or error message
    """
    if not search_dir:
        return ('ok', None)
    dirs = os.walk(search_dir)                                  # liste l'arborescence
    filepath = None
    if os.path.isfile(search_dir):                              # search_dir IS the file
        filepath = search_dir
    elif os.path.isdir(search_dir):                             # search_dir exists
        for dir in dirs:
            if file_name in dir[2]:
                filepath = os.path.join(dir[0], file_name)
                break
    else:                                                       # search_dir doesn't exists
        return (False,  f"directory {search_dir!r} not found.")
    if not filepath:                                            # file not found
        return (False,  f"file {file!r} not found.")
    return (True, filepath)


class Col:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def __usage():
    """
    Help function with argument parser.
    """
    description = ''.join((
                        " =================================\n",
                        " Tool for Reindeer Web App.\n",
                        " Input:\n",
                        "   - Mandatory: fof_unitigs.txt file built to launch Reindeer --index\n",
                        "   - Optional: multiqc_data.json, multiqc_fastq.txt",
                        " or multiqc_general_stats.txt from multiqc output.\n",
                        "               This add in results values to normalize Reindeer query counts.\n"
                        " Output:\n",
                        "   - stdin, or file with '-o' option with Ordered Sample list.\n",
                        " =================================\n",
                        ))

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("fof",
                        help="fof_unitigs.txt file.",
                        metavar=('fof_unitigs.txt'),
    )
    # ~ parser.add_argument("bcalm_log_dir",
                        # ~ help="bcalm log directory (if not used, the number of kmers found will not be recorded and the normalization will not be applied).",
                        # ~ metavar="bcalm/log/directory",
                        # ~ nargs='?',
                        # ~ default=None,
    # ~ )
    parser.add_argument("multiqc_dir",
                        help="Directory of multiqc results.",
                        metavar="multiqc_dir",
                        nargs='?',
                        default=None,
    )
    parser.add_argument('-k', '--kmer-len',
                        type=int,
                        help="kmer length (default: 31)",
                        default=31,
    )
    '''
    parser.add_argument('-s', '--suffix',
                        help="suffix of bcalm log file (default: _bcalm.log)",
                        default='_bcalm.log',
    )
    parser.add_argument('-p', '--prefix',
                        help="prefix of bcalm log file (no prefix by default)",
                        default='',
    )
    '''
    parser.add_argument('-o', '--output',
                        help="output file (default: stdin)",
    )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
    )
    ### Go to "usage()" without arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    __main()
