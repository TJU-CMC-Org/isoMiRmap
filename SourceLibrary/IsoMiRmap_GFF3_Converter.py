#!/usr/bin/env python3
"""IsoMiRmap_GFF3_Converter.py: Converts files from IsoMiRmap table format into GFF3 format. Requires Python 3.5 or later"""

__version__ = "1.0.1"

import sys
import argparse

__author__ = "Jeffrey Ma"
__copyright__ = "Copyright 2020, Thomas Jefferson University - Computational Medicine Center"
__credits__ = ["Jeffrey Ma", "Phillipe Loher"]
# __license__ = "TBD"
__maintainer__ = "Phillipe Loher"
__email__ = "phillipe.loher@jefferson.edu"
__status__ = "Production"

# Checks the version of python being used as the first thing this script does
if sys.version_info[0] < 3 or (sys.version_info[0] >= 3 and sys.version_info[1] < 5):
    sys.stderr.write("Exception: Must be using Python 3.5 or later\n")
    sys.exit(1)


def check_header(check_file):
    """
    Checks the header of the input file
    :param check_file: The input file being checked
    """
    next(check_file)
    next(check_file)
    next(check_file)
    next(check_file)
    next(check_file)
    next(check_file)
    line = check_file.readline()
    if "License Plate" not in line:
        print("Exception: " + check_file.name + " is not in valid format. Missing the header line", file=sys.stderr)
        sys.exit(1)


def header(name, file_name, arguments, args):
    """
    Writes the header required for the GFF3 standard
    :param name: The alias version being used
    :param file_name: The filename of the items being converted. Used to get COLDATA name
    """
    if type(name) == list:
        name = name[0]

    file_name = file_name.split('.')[0]

    print("## mirGFF3. VERSION 1.1");
    print("## Tool runtime arguments: %s", (arguments));
    print("## source-ontology: " + args.mt_name + ' (' + args.mt_rev + ') - ' + args.mt_description);
    print("## COLDATA: " + file_name)
    print("## Filter: PASS - Is ok")
    print("## Filter: EXCLUSIVE - Sequence that exists only within our miR-space (1 or more times) "
          "and not elsewhere in the assembly")
    print("## Filter: AMBIGUOUS - Sequence that exists within our miR-space (1 or more times) "
          "and elsewhere in the assembly")
    print("## Filter: WITHIN_REPEATMASKERISLAND - Which island the sequence is also contained within")


def print_output(output):
    """
    Prints the 3D list of outputs to STDOUT
    :param output: The output list
    """
    for entry in output:
        for column in entry:
            print(column, end='')
        print(';')


def start_conversion(read_file, name, mt_name):
    """
    Begin the conversion process
    :param read_file: The file being read
    :param name: Which alias version to use
    """
    final_output = []
    with read_file as f:
        for line in f:
            # Strip off any trailing characters
            line = line.rstrip()
            if "Ambiguous" in f.name or "ambiguous" in f.name:
                expanded_list = pre_process(line, True)
            else:
                expanded_list = pre_process(line, False)
            for item in expanded_list:
                final_output.append(process(item, name, mt_name))

    print_output(final_output)


def find_index_partial(input_string, search_list):
    """
    Finds the existence of a meta data in the meta data list containing the string
    :param input_string: The string to be searched for
    :param search_list: The list of meta data being searched
    :return: Index in list that the item exists at
    """
    index = 0
    for item in search_list:
        if input_string in item[:2]:
            return index
        else:
            index += 1
    return -1


def pre_process(line, ambiguous):
    """
    This pre-processor splits all the hairpins up into unique entries in a 3D list along with their associated variants
    :param line: Each line of the file
    :param ambiguous: Is the file the ambiguous file or not
    """
    token = line.split('\t')
    if len(token) < 10:
        # There's only 8 columns, so it's missing the island info
        if ambiguous:
            token.append('ambiguous')
        else:
            token.append('exclusive')
    else:
        if ambiguous:
            token[9] = 'ambiguous,' + token[9]
        else:
            token[9] = 'exclusive,' + token[9]

    col7 = token[7].split(', ')
    token.append(str(len(col7)))

    variant = token[8].split('], ')

    # Remove trailing ] on the last entry
    variant[len(variant) - 1] = variant[len(variant) - 1][:-1]

    all_entries = []

    for precursors in range(len(col7)):
        current_entry = []
        for i in range(11):
            if i == 7:
                # Precursor
                current_entry.append(col7[precursors])
            elif i == 8:
                # Meta-data
                if variant[precursors] == '[]' or variant[precursors] == '[':
                    # No meta data
                    current_entry.append('NA')
                else:
                    # The [1:] chops off the '[' part of the meta data
                    current_entry.append(variant[precursors][1:])
            else:
                current_entry.append(token[i])
        all_entries.append(current_entry)

    return all_entries


def variant_calc(five_prime, three_prime, offset3p):
    """
    Calculates the variant information based off the meta data offsets
    :param five_prime: The five prime offset
    :param three_prime: The three prime offset
    :return: The string to be printed
    """
    variant_list = [];
    if five_prime != '0':
        variant_list.append ('iso_5p:' + five_prime);
    if three_prime != '0':
        variant_list.append ('iso_3p:' + three_prime);
    if offset3p > 0:
       variant_list.append ('iso_add3p:' + str (offset3p));
    if len (variant_list) == 0:
        return 'NA'

    return ','.join (variant_list);


def process(entry, name, mt_name):
    """
    Stores all output into a list to be printed later
    :param entry: The entry from the list that contains all the information
    :param name: Which alias version to use
    """
    output = []
    col7 = entry[7].split('&')
    # Col1 seqID
    output.append(col7[0] + '\t')

    # Col2 Source
    output.append(mt_name + "\t")

    # Col3 Type
    output.append("ref_miRNA\t")

    """
    Col4 and 5 Start/End position on chromosome
    Splits col7's third item by their periods. For reference, this is what we're working with at the start:
    8|-|140732581|140732649@7.27.21
    """
    seq_info = col7[2].split('.')

    # The first item is the offset, but needs to be split from the sequence index
    seq_index_info = seq_info[0].split('@')
    strand_info = seq_index_info[0].split('|')

    # Change col 1 to chromosome
    output[0] = strand_info[0] + '\t'

    strand = strand_info[1]
    offset = int(seq_index_info[1])

    # The third item is the length
    string3p = seq_info[2].split ('(');
    length = int(string3p[0])
    offset3p = 0;
    if (len (string3p) > 1):
       offset3p = int (string3p[1][1]); # e.g. only get the 1 in +1U)

    if strand == '+':
        start = int(strand_info[2])
        start_index = start + (offset - 1)
        end_index = start_index + (length - 1) + offset3p  # Subtract once since it's inclusive
    else:
        start = int(strand_info[3])
        start_index = start - (offset - 1) - (length - 1);  # Subtract 2 since it's two inclusive counts
        end_index = start_index + (length - 1)
        start_index -= offset3p;

    output.append(str(start_index) + '\t' + str(end_index) + '\t')

    # Col6 Score
    output.append('0\t')

    # Col7 Positive strand or negative strand
    output.append(strand + '\t')

    # Col8 Phase. Not applicable to us
    output.append('.\t')

    # Col9 Attributes
    output.append('Read=' + entry[1] + '; ')

    # RPM selection dependent
    # output.append(entry[2] + '; ')

    output.append('id=' + entry[0] + '; ')
    output.append('Name=NA; ')
    output.append('Parent=' + col7[0] + '; ')
    output.append('Variant=NA; ')

    if entry[8] != 'NA':
        meta_data = entry[8].split(', ')

        if name == 'miRCarta':
            index = find_index_partial('m-', meta_data)

            if index > -1:
                raw_meta = meta_data[index]
                # Chop off anything after the second &
                # output[output.index('Name=NA; ')] = 'Name=' + '&'.join(raw_meta.split('&', 2)[:2]) + '; '

                # Chop off anything after the first &
                output[output.index('Name=NA; ')] = 'Name=' + raw_meta.split('&', 1)[0] + '; '

                variant = variant_calc(raw_meta.split('|')[4], raw_meta.split('|')[5].split('(')[0], offset3p)
                variant_index = output.index('Variant=NA; ')
                output[variant_index] = 'Variant=' + variant + '; '

                if variant != 'NA':
                    output[2] = 'isomiR\t'
            else:
                output[2] = 'pre_miRNA\t'

        elif name == 'TJU':
            index = find_index_partial('TJ', meta_data)

            if index > -1:
                raw_meta = meta_data[index]
                output[output.index('Name=NA; ')] = 'Name=' + raw_meta.split('&', 1)[0] + '; '

                variant = variant_calc(raw_meta.split('|')[1], raw_meta.split('|')[2].split('(')[0], offset3p)
                variant_index = output.index('Variant=NA; ')
                output[variant_index] = 'Variant=' + variant + '; '

                if variant != 'NA':
                    output[2] = 'isomiR\t'
            else:
                output[2] = 'pre_miRNA\t'

        else:
            # Only one left is miRBase
            index = find_index_partial('MI', meta_data)

            if index > -1:
                raw_meta = meta_data[index]
                insert_index = output.index('Name=NA; ')
                output[insert_index] = 'Name=' + raw_meta.split('&')[1] + '; '

                output.insert(insert_index + 1, 'Alias=' + '&'.join(raw_meta.split('&', 2)[:2]) + '; ')

                variant = variant_calc(raw_meta.split('|')[1], raw_meta.split('|')[2].split('(')[0], offset3p)
                variant_index = output.index('Variant=NA; ')
                output[variant_index] = 'Variant=' + variant + '; '

                if variant != 'NA':
                    output[2] = 'isomiR\t'
            else:
                output[2] = 'pre_miRNA\t'
    else:
       output[2] = 'pre_miRNA\t'

    cigarInsertion = "";
    if offset3p > 0:
       cigarInsertion = str (offset3p) + 'I';
       
    output.append('Cigar=' + str(length) + 'M' + cigarInsertion + '; ')

    output.append('Hits=' + entry[10] + '; ')

    output.append('Genomic=' + strand_info[0] + ':' + str(start_index) + '-' + str(end_index) + '(' + strand + '); ')

    output.append('Expression=' + entry[3] + '; ')

    output.append('normalized=RPM_sameType:' + entry[4] + ',RPM_inputFile:' + entry[5] + ',RPM_custom:'
                  + entry[6] + '; ')

    # Filters
    filters = entry[9].split(',', 1)
    if 'ambiguous' in filters:
        output.append('Filter=AMBIGUOUS')
    else:
        output.append('Filter=PASS,EXCLUSIVE')

    if len(filters) > 1:
        island_filters = filters[1].split(',')
        for island in island_filters:
            output.append(',WITHIN_REPEATMASKERISLAND:' + island)

    return output


def arg_process(arguments):
    parser = argparse.ArgumentParser(
        description='Converts files from IsoMiRmap table format into GFF3 format. Requires Python 3.5 or later.',
        allow_abbrev=True)

    parser.add_argument('exclusive', type=argparse.FileType('r'), help='Exclusive IsoMiRmap table');
    parser.add_argument('ambiguous', type=argparse.FileType('r'), help='Ambiguous IsoMiRmap table');
    parser.add_argument('--mt_name', '--n', type=str, 
                        help='MappingTable bundle name used during mapping (MAPPINGTABLES_NAME).', required=True)
    parser.add_argument('--mt_rev', '--r', type=str, 
                        help='MappingTable bundle revision used during mapping (MAPPINGTABLES_REV).', required=True)
    parser.add_argument('--mt_description', '--d', type=str, 
                        help='MappingTable bundle description used during mapping (MAPPINGTABLES_DESCRIPTION).', required=True)
    parser.add_argument('--meta', '--m', type=str, choices=['miRCarta', 'TJU', 'miRBase'], default='miRCarta',
                        help='Which meta data entry type to use for the output for the name parameter. '
                             'Default is miRCarta. Case sensitive.', required=True);

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(arguments)
    conversion(args, arguments)


def conversion(args, arguments):
    exclusive = args.exclusive
    ambiguous = args.ambiguous
    name = args.meta

    # Error checks to make sure the header is present
    check_header(exclusive)
    check_header(ambiguous)

    # Error checks to make sure the prefixes on both files match
    exclusive_check = str(exclusive.name).split('-')
    ambiguous_check = str(ambiguous.name).split('-')

    exclusive_check.remove('exclusive')
    ambiguous_check.remove('ambiguous')

    if exclusive_check != ambiguous_check:
        for i in range(len(exclusive_check)):
            if exclusive_check[i] != ambiguous_check[i]:
                print('Exception: Unexpected file names.', file=sys.stderr, end=' ')
                print("Exclusive's file name contains " + exclusive_check[i] +
                      " but ambiguous' file name contains " + ambiguous_check[i],
                      file=sys.stderr)
                sys.exit(1)

    # Begin conversion
    header(name, str(exclusive.name), arguments, args)
    start_conversion(exclusive, name, args.mt_name)
    start_conversion(ambiguous, name, args.mt_name)

    # Making sure everything is closed
    exclusive.close()
    ambiguous.close()


# Initiates the program
if __name__ == "__main__":
    arg_process(sys.argv[1:])
