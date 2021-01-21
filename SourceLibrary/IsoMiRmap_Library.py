#!/usr/bin/env python3
"""IsoMiRmap.py: Use this program to generate isomiR profiles. Requires Python 3.5 or later"""

__version__ = '1.0.0'

import sys
import argparse
import os
import io
import hashlib
import re
from collections import OrderedDict
from pathlib import Path;
import datetime
import re
import SourceLibrary.MINTplates as MINT
import SourceLibrary.IsoMiRmap_GFF3_Converter as Converter
import gzip

__info__= 'Thomas Jefferson University - Computational Medicine Center'
__credits__ = ['Phillipe Loher', 'Jeffrey Ma']
__maintainer__ = 'Phillipe Loher'
__email__ = 'phillipe.loher@jefferson.edu'
__status__ = 'Development'

# globals
scriptversion = 'IsoMiRmap_v5' # change version if output changes in anyway (new logic, new functionality, or even if the format of the output changes)
DEFAULT_OUTPREFIX = 'output'
DEFAULT_FRAGMENTTYPE = 'isomiR'
LICENSEPLATEPREFIX = 'iso'
DEFAULT_MAPPINGTABLES="MappingBundles/miRCarta";
stat_totalstartingreads = 0
min_fragment_len = None;
mappingtables_name = None
mappingtables_rev = None
mappingtables_description = None

# MD5sum and table file names, tables must match this to make sure they are paired correctly
fn_lookuptable = None
md5sum_lookuptable = "na"
fn_splicedsequences = None
md5sum_splicedsequences = "na"
fn_otherannotations = None
md5sum_otherannotations = "na"
fn_metacoordinates = None
md5sum_metacoordinates = "na"
fn_snptable = None
md5sum_snptable = "na"

# regex compiles
ALPHANUMERIC = re.compile('^[&-.|+_a-z0-9]+$', re.IGNORECASE)
VALID = re.compile('^[ATCGN]+$', re.IGNORECASE)

# open a file for reading, handle if it has gz extension
def open_read_gz (arg_fn, MD5_check):
    retfh = None
    if '.gz' in str (arg_fn):
        try:
            retfh = io.BufferedReader(gzip.open(str (arg_fn), 'rb'))
        except IOError:
            print('Error, exiting: file ' + str (arg_fn) + ' could not be opened', file=sys.stderr)
            sys.exit(1)
    else:
        try:
            retfh = open(str (arg_fn), 'rb')
        except IOError:
            print('Error, exiting: file ' + str (arg_fn) + ' could not be opened', file=sys.stderr)
            sys.exit(1)
           
    if MD5_check != None: 
       md5check = hashlib.md5(retfh.read()).hexdigest();
       if MD5_check != md5check:
           print('Error, exiting: md5sum defined in tables.cfg does not match for {}: {} vs {}'.format (str (arg_fn), MD5_check, md5check), file=sys.stderr);
           sys.exit(1)
       retfh.seek (0)
    return retfh;


def process_args():
    """
    Process arguments
    :return: processed arguments
    """

    # Need the path of the script because files are being taken in as strings and not directly opened.
    # Path is appended to the string name to get actual path to file
    path = os.path.dirname(os.path.realpath(sys.argv[0]))

    if path == os.getcwd():
        path = ''
    else:
        if os.name == 'posix':
            path += '/'
        elif os.name == 'nt':
            path += '\\'

    from argparse import RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description='Use this program to generate isomiR profiles. Requires Python 3.5 or later.\n'
                    'OUTPUT:\n'
                    '- Plain text and HTML file pairs for exclusive isomiRs '
                    '(*.exclusive-isomiRs.expression.*). '
                    'IsomiRs that exist exclusively within miR-space (higher confidence) will be reported in this file.  '
                    'Also, their counterparts with 3p additions (e.g. uridylation) will be included.  '
                    'RPM and annotation information included.\n'
                    '- Plain text and HTML file pairs for ambiguous isomiRs '
                    '(*.ambiguous-isomiRs.expression.*). Same as above but for isomiRs that exist both within miR-space '
                    'and also elsewhere in the genome.\n'
                    '- Variant-containing isomiRs '
                    '(*.snps-isomiRs.expression.*)\n'
                    '- High level mapping stats are also generated seperately for exclusive, non-exclusive, and variant containing isomiRs '
                    '(*.countsmeta.txt)\n'
                    '- Output in mirGFF3 format '
                    '(*.gff3)', 
        allow_abbrev=False, formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('trimmedfastqfile', type=str,
                        help='This is a required argument. The file “trimmedfastqfile” contains the sequencer reads'
                             ' from the Short-RNA dataset. Any trimming (e.g. quality and adapter trimming)'
                             ' must be done prior to running this tool. If “trimmedfastqfile” ends in “.gz” it will be'
                             ' treated as a gzipped FASTQ file. The FASTQ file contains four lines per read (more info'
                             ' here: https://en.wikipedia.org/wiki/FASTQ_format). Color-space reads are not supported.');
    parser.add_argument('--m', '--mappingbundle', type=str, default=DEFAULT_MAPPINGTABLES, dest='mappingtables',
                        help='If not specified, “mappingbundle” will be set to "' + DEFAULT_MAPPINGTABLES + '". This is the'
                             ' relative or absolute path to the mapping bundle that will be used.  Out of the box, we'
                             ' provide mapping bundles for both miRCarta and miRBase for Homo Sapiens.');
    parser.add_argument('--p', '--outputprefix', type=str, default=DEFAULT_OUTPREFIX, dest='outputprefix',
                        help='If not specified, “outputprefix” will be set to "' + DEFAULT_OUTPREFIX + '" and all'
                             ' output files will be generated in the current working directory. The “outputprefix” is a'
                             ' string that will be prepended in the output files that ' + scriptversion + ' will generate. The'
                             ' string is meant to serve as a mnemonic for the user. “outputprefix” can optionally include'
                             ' relative or absolute directory paths (if the directories exist) to save the output to a'
                             ' different directory.')
    parser.add_argument('--d', '--customRPM', type=int, default=None, dest='customrpm',
                        help='This is an optional argument. The value “customRPM” is meant to be used as alternative'
                             ' denominator when computing the RPM abundance of an isomiR. When this parameter is'
                             ' defined by the user, an additional column in the output files will be populated with'
                             ' RPM values that have been calculated using this value in the denominator – i.e. these'
                             ' values will be equal to raw reads/<customRPM>*1,000,000. A common value to use here is'
                             ' the original number of sequenced reads prior to quality- and adaptor-trimming.')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if not os.path.exists(str (Path (args.mappingtables) / "tables.cfg")):
        print('Error, exiting: Config file tables.cfg not found in directory: ' + args.mappingtables, file=sys.stderr)
        sys.exit(1)
        
    return args

def load_configtable (fn_configtable):
    """
    Loads lookup table entries into their respective dictionaries
    :param fn_configtable: The name of the config table file
    """
    print('Loading Config Table')
    global mappingtables_name;
    global mappingtables_rev;
    global mappingtables_description;
    global md5sum_lookuptable;
    global fn_lookuptable;
    global md5sum_splicedsequences;
    global fn_splicedsequences;
    global md5sum_otherannotations;
    global fn_otherannotations;
    global md5sum_metacoordinates;
    global fn_metacoordinates;
    global md5sum_snptable;
    global fn_snptable;

    try:
        open_lookup = open(str (fn_configtable), 'r')
    except IOError:
        print('Error, exiting: ' + str (fn_configtable) + ' not found', file=sys.stderr)
        sys.exit(1)

    for line in open_lookup:
        read_line = line.split ("#")[0].strip (); # remove any comments from line if they exist
        split_value = read_line.split ("=");
        if len (split_value) > 1:
           thekey = split_value[0].strip ().upper ();
           thevalue = split_value[1].strip ();
           md5val = "na"
           md5split = thevalue.split ("MD5SUM");
           if len (md5split) > 1:
              md5val = split_value[2].strip ();
              thevalue = md5split[0].strip ();
    
           if thekey == "MAPPINGTABLES_NAME":
              mappingtables_name = thevalue;
           elif thekey == "MAPPINGTABLES_REV":
              mappingtables_rev = thevalue;
           elif thekey == "MAPPINGTABLES_DESCRIPTION":
              mappingtables_description = thevalue;
           elif thekey == "FRAGMENTS":
              md5sum_lookuptable = md5val
              fn_lookuptable = thevalue
           elif thekey == "HAIRPINSEQUENCES":
              md5sum_splicedsequences = md5val
              fn_splicedsequences = thevalue
           elif thekey == "OTHERANNOTATIONS":
              md5sum_otherannotations = md5val
              fn_otherannotations = thevalue
           elif thekey == "METACOORDINATES":
              md5sum_metacoordinates = md5val
              fn_metacoordinates = thevalue
           elif thekey == "SNPS":
              md5sum_snptable = md5val
              fn_snptable = thevalue
           
    if (fn_lookuptable == None):
        print('Error, exiting: FRAGMENTS not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
    if (fn_splicedsequences == None):
        print('Error, exiting: HAIRPINSEQUENCES not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
    if (fn_metacoordinates == None):
        print('Error, exiting: METACOORDINATES not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
        
    if (mappingtables_name == None):
        print('Error, exiting: MAPPINGTABLES_NAME not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
    if (mappingtables_rev == None):
        print('Error, exiting: MAPPINGTABLES_REV not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
    if (mappingtables_description == None):
        print('Error, exiting: MAPPINGTABLES_DESCRIPTION not set in tables.cfg', file=sys.stderr);
        sys.exit(1)
       
    open_lookup.close()

def load_splicedfasta(sequence):
    """
    Loads spliced sequences into a dictionary
    :param sequence: The path/name of the sequence file
    :return: The fully loaded disctionary
    """
    print('Loading Hairpin sequences')
    hash_splicedseqs = {}

    open_sequence = open_read_gz (sequence, md5sum_splicedsequences);
    header = 1
    hairpinname = None
    for line in open_sequence:
        if header:
            # Strip off any \n or trailing whitespaces
            header = line.decode().rstrip()
            if header[0] != '>':
                print('Error, exiting: Expected FASTA header but received ' + header + ' instead.', file=sys.stderr)
                sys.exit(1)

            # Remove the > from the header name to get the hairpin name and assign to variable
            hairpinname = header[1:]

            if not ALPHANUMERIC.match(hairpinname):
                print('Error, Only characters [&-.|+_A-Za-z0-9] allowed in FASTA label but received ' + header +
                      '. Note, no whitespace is allowed', file=sys.stderr)
                sys.exit(1)
            header = 0
        else:
            sequence = line.decode().rstrip()
            if not VALID.match(sequence):
                print('Error, exiting: Only characters [ATGCN] allowed in FASTA sequence but received '
                      + sequence + ' instead', file=sys.stderr)
                sys.exit(1)

            # Add to hash
            hash_splicedseqs[hairpinname] = sequence
            header = 1

    open_sequence.close()

    return hash_splicedseqs


def load_otherannotations(annotations):
    """
    Loads annotations into dictionary
    :param annotations: The annotation file name
    :return: Loaded annotations dictionary
    """
    print('Loading OtherAnnotations')
    hash_otherannotations = {}

    open_annotations = open_read_gz (annotations, md5sum_otherannotations);
    for line in open_annotations:
        # Strip off any \n or trailing whitespaces
        read_line = line.decode().rstrip()
        split_lines = read_line.split('\t')

        # Add to hash
        if len(split_lines) == 1:
            linevalue = ''
        else:
            linevalue = split_lines[1]
        hash_otherannotations[split_lines[0]] = linevalue

    open_annotations.close()
    return hash_otherannotations
 
def load_snps(annotations):
    """
    Loads annotations into dictionary
    :param annotations: The annotation file name
    :return: Loaded annotations dictionary
    """
    print('Loading SNP Table')
    hash_snps = {}

    open_snps = open_read_gz (annotations, md5sum_snptable);
    for line in open_snps:
        # Strip off any \n or trailing whitespaces
        read_line = line.decode().rstrip()
        split_lines = read_line.split('\t')

        # Add to hash
        hash_snps[split_lines[0]] = split_lines[1:];

    open_snps.close()
    return hash_snps;


def load_metacoordinates(metacoordinates):
    """
    Loads meta coordinates into dictionary
    :param metacoordinates: Meta coordinate file name
    :return: Loaded meta coordinates dictionary
    """
    print('Loading Metadata coordinates')
    hash_metacoordinates = {}

    open_metacoordinates = open_read_gz (metacoordinates, md5sum_metacoordinates);
    for line in open_metacoordinates:
        # Strip off any \n or trailing whitespaces
        read_line = line.decode().rstrip()
        split_lines = read_line.split('\t')

        # Add to hash. Key is contig + strand
        hashkey = split_lines[0] + split_lines[1]
        if hashkey not in hash_metacoordinates:
            hash_metacoordinates[hashkey] = []
        savearray = [split_lines[2], split_lines[3], split_lines[4]]
        hash_metacoordinates[hashkey].append(savearray)

    open_metacoordinates.close()
    return hash_metacoordinates


def load_lookup_table(lookuptable):
    """
    Loads lookup table entries into their respective dictionaries
    :param lookuptable: The name of the lookup table file
    :return: A list of the two dictionaries
    """
    print('Loading Lookup Table')
    global min_fragment_len
    hash_exclusive = {}
    hash_notexclusive = {}

    open_lookup = open_read_gz (lookuptable, md5sum_lookuptable);
    for line in open_lookup:
        # Strip off any \n or trailing whitespaces
        read_line = line.decode().rstrip()
        line_split = read_line.split('\t')
        
        if min_fragment_len == None or len (line_split[0]) < min_fragment_len:
           min_fragment_len = len (line_split[0]); # keep track of the minimum fragment size inside the lookup tasble

        # Note: a Y in col 2 means the fragment sequence is exclusive to miRNA-space and a N means it's not exclusive
        if line_split[1] == 'Y':
            hash_exclusive[line_split[0]] = 0
        elif line_split[1] == 'N':
            hash_notexclusive[line_split[0]] = 0
        else:
            print('Error, exiting: Lookup value of ' + line_split[1] + ' in ' + lookuptable + ' is not recognized',
                  file=sys.stderr)

    open_lookup.close()
    return [hash_exclusive, hash_notexclusive]


def add_annotations(array_addto, isomirfull, isomircut, hash_hairpins):
    """
    Creates the hairpin annotation information for
    :param isomirseq:
    :param hash_hairpins:
    :return:
    """
    array_annotationoutput = []
    onematch = False; # for error checking, make sure isomircut is at least found once or else something is wrong
    isomircut_len = len(isomircut)
    largestmod = 0;

    # Look for the isomiR seq across all hairpinspliced sequences
    for hairpinname in hash_hairpins.keys():
        hairpinseq = hash_hairpins[hairpinname]
        offset = 0
        result = hairpinseq.find(isomircut, offset)
        while result != -1:
            onematch = True;
            matchprev = hairpinname + '@' + str(result + 1) + '.'; # this is a match key to see if it was found in the same hairpin at the same start location
            
            addme=True; # add the hairpin location if it hasn't been added already (potentially between T-trimming)
            for theannotation in array_addto:
               if theannotation.startswith(matchprev):
                  addme=False;
                  break;
            
            if addme == True: # add outside of the iterator above since adding while searching may cause issues
               mod3p = '';
               if (len (isomirfull) > isomircut_len):
                  modsize = len(isomirfull) - isomircut_len;
                  lastnt = isomirfull[-1:];
                  if (lastnt == 'T'):
                    lastnt = 'U';
                  mod3p = '(+' + str (modsize) + lastnt + ')'
                  if (modsize > largestmod):
                     largestmod = modsize;
               newadd = hairpinname + '@' + str(result + 1) + '.' + str(result + isomircut_len) + '.' + str(isomircut_len) + mod3p
               array_addto.append (newadd)
               
            offset = result + 1
            result = hairpinseq.find(isomircut, offset)

    if not onematch:
        print("Error, exiting: no annotations in sequence file " + fn_splicedsequences + " found for isomiR "
              + isomircut, file=sys.stderr)
        sys.exit(1)
        
    return largestmod



def get_meta_intersects(metaannotations_hash, annotation_list):
    """
    Get intersections of meta coordinates
    :param metaannotations_hash: The dictionary containing all the meta annotations
    :param annotation_list: A list of annotation entries, such as hsa-1-52.1&WithFlank&17|+|59841267|59841338@7.28.22
    :return: An array of arrays containing all meta coordinate intersections
    """
    retaa = []

    for annotationentry in annotation_list:
        ret = []
        # example of annotationentry: hsa-1-52.1&WithFlank&17|+|59841267|59841338@7.28.22
        amp_split = annotationentry.split('&')
        locations_split = amp_split[2].split('@')
        anno_coords = locations_split[0].split('|')
        anno_offsets = list(map(int, locations_split[1].split('(')[0].split('.')))
        
        postmod = ''; # lets carry over the hairpin post-modification annotations to any mature intersect annotations
        if (len (locations_split[1].split('(')) == 2):
           postmod = '(' + locations_split[1].split('(')[1];
           
        isomir_coord_start = -1
        isomir_coord_end = -1

        # Converts the string numbers into actual ints
        anno_coords = [int(x) if x.isnumeric() else x for x in anno_coords]

        if anno_coords[1] == '+':
            isomir_coord_start = anno_coords[2] + anno_offsets[0] - 1
            isomir_coord_end = anno_coords[2] + anno_offsets[1] - 1
        elif anno_coords[1] == '-':
            isomir_coord_start = anno_coords[3] - anno_offsets[1] + 1
            isomir_coord_end = anno_coords[3] - anno_offsets[0] + 1
        else:
            print("Exiting... invalid isomiR strand " + anno_coords[1])
            exit(1)

        # only search for coordinates with the same coordinate/strand to speed things up.
        # This could be made much faster as well with pre-sorting during loading etc. (future)
        keylookup = str(anno_coords[0]) + anno_coords[1]
        if keylookup in metaannotations_hash:
            # if coordinates exist for the same annotation location, lets see if there's intersection
            for coordinate_intersects in metaannotations_hash[keylookup]:

                # This commented out section was preserved from the initial conversion from perl
                # print("Comparing {} with {}".format('-'.join(anno_coords), '-'.join(coordinate_intersects)))
                # print("coordinates {} and {}".format(coordinate_intersects[0], coordinate_intersects[1]))
                intersect_start = int(coordinate_intersects[0])
                intersect_end = int(coordinate_intersects[1])

                # report the coordinate if there's any intersection
                if intersect_start <= isomir_coord_start <= intersect_end \
                        or intersect_start <= isomir_coord_end <= intersect_end \
                        or (isomir_coord_start < intersect_start and isomir_coord_end > intersect_end):
                    # now get isomiR offset from coordinates (note: this calculation is dependent on strand).
                    # Upstream offsets should be negative
                    offset5p = ''
                    offset3p = ''
                    if anno_coords[1] == '+':
                        offset5p = isomir_coord_start - intersect_start
                        offset3p = isomir_coord_end - intersect_end
                    elif anno_coords[1] == '-':
                        offset5p = intersect_end - isomir_coord_end
                        offset3p = intersect_start - isomir_coord_start

                    # To be consistent with our previous naming convention, lets add the plus sign for positive numbers.
                    # Negative numbers will already contain the negative
                    if offset5p > 0:
                        offset5p = '+' + str(offset5p)
                    if offset3p > 0:
                        offset3p = '+' + str(offset3p)
                    record_result = coordinate_intersects[2] + '&offsets|' + str(offset5p) + '|' + str(offset3p) + postmod
                    ret.append(record_result)
        if (ret == []):
           retaa.append(['No mature intersections'])
        else:
           retaa.append(ret)
         
    return retaa


def write_argsheader(file_towrite, html):
    """
    Writes the program arguments to file as first line
    :param file_towrite: File to write to
    :param html: If the file being written to is a html file
    """
    today = datetime.datetime.now()

    if not html:
        file_towrite.write('## ')

    file_towrite.write(scriptversion + '\t'
                       + str(today.month) + '/' + str(today.day) + '/' + str(today.year) + '\t')
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i]
        if '\\' in arg:
            arg.replace('\\', '/')
        if '/' in arg:
            arg = arg.rsplit('/', 1)[1]
        if i == len(sys.argv) - 1:
            if html:
                file_towrite.write(arg + '<br /><br />' + '\n')
            else:
                file_towrite.write(arg + '\n')
        else:
            file_towrite.write(arg + '\t')


def create_output(read_hash, annotation_hash, otherannotations_hash, largest_postmod,
                  total_frags_in_file, seqcategory, metaannotations_hash, output_prefix, customrpm):
    """
    Actually generates the output file
    :param read_hash: Dictionary of all the read sequences
    :param annotation_hash: Dictionary of all annotations
    :param otherannotations_hash: Dictionary of all otherannotations
    :param largest_postmod: Dictionary of the largest postmod stretch for a given sequence
    :param total_frags_in_file: How many fragments were in the read file
    :param seqcategory: Sequence category of read fragment, either exclusive or ambiguous
    :param metaannotations_hash: Dictionary of all metaannotations
    :param output_prefix: Prefix to be used for the final file. User can also specify the path along with the prefix
    :param customrpm: CustomRPM to be used for division. Column is filled with 'na' if not specified
    :return: The filename of the text file created. For use in GFF3 conversion.
    """
    # output the fragments in decreasing order of expression
    # double counting is not a problem because we are dealing with the raw reads
    ret_file_name = output_prefix + '-' + scriptversion + '-' + seqcategory + '.expression.txt'
    print("Creating output file: " + ret_file_name)

    expression_file = open(ret_file_name, 'w+')
    expression_file.write ('## Table of ' + seqcategory + '.\n## ' + mappingtables_name + ' (' + mappingtables_rev + '): ' + mappingtables_description + '\n')
    write_argsheader(expression_file, False)
    expression_file.write ("## RPM* - Using [{}] as denominator - total number of reads in this output file.\n"
                           "## RPM** - Using [{}] as denominator - total number of reads in FASTQ file.\n"
                           "## RPM*** - Using [{}] as denominator - custom denominator passed in with -d parameter.\n".format(total_frags_in_file,
                                                                                                                             stat_totalstartingreads,
                                                                                                                             customrpm if customrpm is not None else 'na'))

    expression_file.write("License Plate\tIsomiR sequence\tType\tUnnormalized read counts\tRPM*\tRPM**\tRPM***\t"
                          "Hairpin locations (comma delimited)\tMature"
                          " meta-data (bracket delimited per hairpin)\tRepeatMaskerClassIsland where fragment is fully"
                          " contained (comma delimited)\n");

    filename = output_prefix + '-' + scriptversion + '-' + seqcategory + '.expression.html'
    print("Creating output file: " + filename)

    html_file = open(filename, 'w+')

    html_file.write("<html><head><title>{} expression</title></head>\n".format(seqcategory))
    html_file.write('<style> tr:nth-of-type(odd) { background-color:#ccc; } body { font-size: 18px; } table'
                    ' { font-size: 16px; } </style>\n')
    html_file.write('<body>'
                    'Table of ' + seqcategory + '.  <br />' + mappingtables_name + ' (' + mappingtables_rev + '): ' + mappingtables_description + '<br />\n')
    write_argsheader(html_file, True)
    html_file.write ("<p style=\"font-size:12px; display:inline\">"
                     "RPM* - Using [{}] as denominator - total number of reads in this output file. <br />"
                     "RPM** - Using [{}] as denominator - total number of reads in FASTQ file. <br />"
                     "RPM*** - Using [{}] as denominator - custom denominator passed in with -d parameter.<br />"
                     "</p>".format(total_frags_in_file,
                            stat_totalstartingreads,
                            customrpm if customrpm is not None else 'na'))

    # This commented out section was preserved from the initial conversion from perl
    """html_file.write('<br />\nCreated by the <a target="_blank" href="http://cm.jefferson.edu">Computational Medicine'
                ' Center</a> at <a target="_blank" href="http://www.jefferson.edu/">Thomas Jefferson University</a> '
                'using the MINTmap tool located <a target="_blank" href="http://cm.jefferson.edu/MINTcodes/">here'
                '</a>.<br />\nPlease cite: Loher, P. <i>et al.</i> MINTmap: fast and exhaustive profiling of nuclear '
                'and mitochondrial tRNA fragments from short RNA-seq data. <i>Sci. Rep.</i> 7, 41184; doi: 10.1038/'
                'srep41184 (2017).')"""

    html_file.write('<table style="width=100%; white-space: nowrap;"><tr><td><b>License Plate</b><br />\n(sequence derived)</td>'
                    '<td><b>IsomiR sequence</b></td>'
                    '<td><center><b>Type</b></center></td>'
                    '<td><center><b>Unnormalized<br />\nread counts</b></center></td>'
                    '<td><center><b>RPM*</b></center></td>'
                    '<td><center><b>RPM**</b></center></td>'
                    '<td><center><b>RPM***</b></center></td>'
                    '<td><b>Hairpin locations<br />\n'
                    '(comma delimited)</b></td><td><b>Mature meta-data<br />\n(bracket delimited per '
                    'hairpin)</b></td><td><b>RepeatMaskerClassIsland where fragment is fully contained<br />\n'
                    '(comma delimited)</b></td></tr>');

    read_hash = OrderedDict(sorted(read_hash.items(), key=lambda x: (-x[1], x[0])))

    for key in read_hash:
        if read_hash[key] != 0:
            # Ensure we don't print things with 0 expression
            unnorm_numreads = read_hash[key]
            rpm_1 = round((unnorm_numreads/total_frags_in_file) * 1000000, 2)
            rpm_2 = round((unnorm_numreads/stat_totalstartingreads) * 1000000, 2)
            if customrpm is None:
                rpm_3 = 'na'
            else:
                rpm_3 = round((unnorm_numreads/customrpm)*1000000, 2)
            annotations = ', '.join(annotation_hash[key])
            annotations_html = ',<br />'.join(annotation_hash[key])

            otherannotations = 'na'
            modsize = 0;
            if fn_otherannotations != None:
                try:
                   modsize = largest_postmod[key];
                   if modsize == 0:
                     otherannotations = otherannotations_hash[key]
                   else:
                     otherannotations = otherannotations_hash[key[:-modsize]] # get repeat masker annotations from the smallest possible templated sequence of the isomiR
                      
                except KeyError:
                    print("Using OtherAnnotations but OtherAnnotation for " + key + " is not found, exiting...",
                          file=sys.stderr)
                    exit(1)

            metaannotations = 'na'
            metaarray = []
            arrayofarrays = get_meta_intersects(metaannotations_hash, annotation_hash[key])
            # this helps deliminate which meta-coordinates coorelate with which precursor in the previous column.
            # The meta-data for each precursor will be separated by " & "
            for i in range(len(arrayofarrays)):
                getarray = arrayofarrays[i]
                metaarray.append(', '.join(getarray))

            metaannotations = ', '.join(list(map(lambda x: "[" + x + "]", metaarray)))
            metaannotations_html = ',<br />'.join(list(map(lambda x: "[" + x + "]", metaarray)))
           
            lastnt = key[-1:];
            if (lastnt == 'T'):
                lastnt = 'U';
            expression_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(MINT.convert(key, True, 'iso'),
                                                                                key,
                                                                                'wildtype' if modsize == 0 else ('wildtype+' + lastnt + '(s)'),
                                                                                unnorm_numreads,
                                                                                rpm_1, rpm_2, rpm_3,
                                                                                annotations,
                                                                                metaannotations,
                                                                                otherannotations))
            html_file.write('<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}'
                            '</td><td>{}</td></tr>\n'.format(MINT.convert(key, True, 'iso'),
                                                             key,
                                                            'wildtype' if modsize == 0 else ('wildtype+' + lastnt + '(s)'),
                                                             unnorm_numreads,
                                                             rpm_1, rpm_2, rpm_3,
                                                             annotations_html,
                                                             metaannotations_html,
                                                             otherannotations))

    html_file.write('</table></body></html>\n')

    expression_file.close()
    html_file.close()

    filename = output_prefix + '-' + scriptversion + '-' + seqcategory + '.countsmeta.txt'
    print('Creating output file ' + filename)

    countsmeta = open(filename, 'w+')
    write_argsheader(countsmeta, False)

    countsmeta.write("Total reads in -f input file\tTotal unnormalized reads in " + seqcategory + "\tPercent\n")
    countsmeta.write("{}\t{}\t{}%".format(stat_totalstartingreads,
                                          total_frags_in_file,
                                          total_frags_in_file/stat_totalstartingreads * 100))

    countsmeta.close()

    return ret_file_name
 
def create_snp_output (read_hash, annotation_hash, total_frags_in_file, seqcategory, output_prefix, customrpm):
    """
    Actually generates the output file
    :param read_hash: Dictionary of all the read sequences and their counts
    :param annotation_hash: Dictionary of all annotations
    :param total_frags_in_file: How many fragments were in the read file
    :param seqcategory: Sequence category of read fragment
    :param output_prefix: Prefix to be used for the final file. User can also specify the path along with the prefix
    :param customrpm: CustomRPM to be used for division. Column is filled with 'na' if not specified
    :return: The filename of the text file created. Potentially for use in GFF3 conversion.
    """
    # output the fragments in decreasing order of expression
    # double counting is not a problem because we are dealing with the raw reads
    ret_file_name = output_prefix + '-' + scriptversion + '-' + seqcategory + '.expression.txt'
    print("Creating output file: " + ret_file_name)

    expression_file = open(ret_file_name, 'w+')
    expression_file.write ('## Table of ' + seqcategory + '.\n## ' + mappingtables_name + ' (' + mappingtables_rev + '): ' + mappingtables_description + '\n')
    write_argsheader(expression_file, False)
    expression_file.write ("## RPM* - Using [{}] as denominator - total number of reads in this output file.\n"
                           "## RPM** - Using [{}] as denominator - total number of reads in FASTQ file.\n"
                           "## RPM*** - Using [{}] as denominator - custom denominator passed in with -d parameter.\n".format(total_frags_in_file,
                                                                                                                             stat_totalstartingreads,
                                                                                                                             customrpm if customrpm is not None else 'na'))

    expression_file.write("License Plate\tIsomiR sequence\tUnnormalized read counts\t"
                          "RPM*\tRPM**\tRPM***\t"
                          "SNP id(s)\tHairpin(s)\tReference sequence(s)\tReference sequence license plate(s)\n")

    filename = output_prefix + '-' + scriptversion + '-' + seqcategory + '.expression.html'
    print("Creating output file: " + filename)

    html_file = open(filename, 'w+')

    html_file.write("<html><head><title>{} expression</title></head>\n".format(seqcategory))
    html_file.write('<style> tr:nth-of-type(odd) { background-color:#ccc; } body { font-size: 18px; } table'
                    ' { font-size: 16px; } </style>\n')
    html_file.write('<body>'
                    'Table of ' + seqcategory + '.  <br />' + mappingtables_name + ' (' + mappingtables_rev + '): ' + mappingtables_description + '<br />\n')

    write_argsheader(html_file, True)
    html_file.write ("<p style=\"font-size:12px; display:inline\">"
                     "RPM* - Using [{}] as denominator - total number of reads in this output file. <br />"
                     "RPM** - Using [{}] as denominator - total number of reads in FASTQ file. <br />"
                     "RPM*** - Using [{}] as denominator - custom denominator passed in with -d parameter.<br />"
                     "</p>".format(total_frags_in_file,
                            stat_totalstartingreads,
                            customrpm if customrpm is not None else 'na'))

    html_file.write('<table style="width=100%; white-space: nowrap;"><tr><td><b>License Plate</b><br />\n(sequence derived)</td><td><b>'
                    'IsomiR sequence</b></td><td><center><b>Unnormalized<br />\nread counts</b></center></td>'
                    '<td><center><b>RPM*</b></center></td>'
                    '<td><center><b>RPM**</b></center></td>'
                    '<td><center><b>RPM***</b></center></td>'
                    '<td><b>SNP id(s)</b></td>'
                    '<td><b>Hairpin(s)</b></td>'
                    '<td><b>Reference sequence(s)</b></td>'
                    '<td><b>Reference sequence license plate(s)</b></td>'
                    '</tr>');
                    

    read_hash = OrderedDict(sorted(read_hash.items(), key=lambda x: (-x[1], x[0])))

    for key in read_hash:
        if read_hash[key] != 0:
            # Ensure we don't print things with 0 expression
            unnorm_numreads = read_hash[key]
            rpm_1 = round((unnorm_numreads/total_frags_in_file) * 1000000, 2)
            rpm_2 = round((unnorm_numreads/stat_totalstartingreads) * 1000000, 2)
            if customrpm is None:
                rpm_3 = 'na'
            else:
                rpm_3 = round((unnorm_numreads/customrpm)*1000000, 2)
                
            splitannotations = annotation_hash[key];
           
            # convert reference sequences to license plates as well 
            license_convert = splitannotations[2];
            re_convert = re.findall (r'[ATCG]+', license_convert)
            for convert in re_convert:
               license_convert = license_convert.replace (convert, MINT.convert(convert, True, 'iso'), 1); # only replace first instance just incase

            expression_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(MINT.convert(key, True, 'iso'),
                                                                                key,
                                                                                unnorm_numreads,
                                                                                rpm_1, rpm_2, rpm_3,
                                                                                splitannotations[0],
                                                                                splitannotations[1], 
                                                                                splitannotations[2],
                                                                                license_convert));
            html_file.write('<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}'
                            '</td><td>{}</td><td>{}</td></tr>\n'.format(MINT.convert(key, True, 'iso'),
                                                             key,
                                                             unnorm_numreads,
                                                             rpm_1, rpm_2, rpm_3,
                                                             splitannotations[0].replace ('], ', '], <br />'),
                                                             splitannotations[1].replace ('], ', '], <br />'), 
                                                             splitannotations[2].replace ('], ', '], <br />'),
                                                             license_convert.replace ('], ', '], <br />')));

    html_file.write('</table></body></html>\n')

    expression_file.close()
    html_file.close()

    filename = output_prefix + '-' + scriptversion + '-' + seqcategory + '.countsmeta.txt'
    print('Creating output file ' + filename)

    countsmeta = open(filename, 'w+')
    write_argsheader(countsmeta, False)

    countsmeta.write("Total reads in -f input file\tTotal unnormalized reads in " + seqcategory + "\tPercent\n")
    countsmeta.write("{}\t{}\t{}%".format(stat_totalstartingreads,
                                          total_frags_in_file,
                                          total_frags_in_file/stat_totalstartingreads * 100))

    countsmeta.close()

    return ret_file_name


def create_gff3_output(outputprefix, args):
    """
    Creates the GFF3 files
    :param outputprefix: The prefix to be used for the output
    :param args: Arguments to be passed into the converter
    :return:
    """
    if '--meta' in args:
        gff3_filename = outputprefix + '-' + scriptversion + '-' + DEFAULT_FRAGMENTTYPE + 's.' + \
                        args[args.index('--meta') + 1] + '.gff3'
        print('Creating GFF3 file ' + gff3_filename)
        gff3 = open(gff3_filename, 'w+')

        # Store standard out
        stdout_restore = sys.stdout
        sys.stdout = gff3

        Converter.arg_process(args)
        gff3.close()

        # Restore standard out
        sys.stdout = stdout_restore
    else:
        print('Error in GFF3 args')
        sys.exit(1)


def start_generation():
    """
    Begins the generation of files
    """
    print(sys.argv[0] + ' (' + scriptversion + ') starting')

    # load arguments passed into script and perform argument error checking
    args = process_args()
    load_configtable (Path (args.mappingtables) / "tables.cfg");

    # load hairpin sequences into memory
    # Dictionary for isomiRs that store all possible annotations for the
    # isomiR in miRNA-space: hash(isomiRseq)->(array of strings)
    fastfrag_hairpin_annotations = {}
    fastfrag_largestPostMod = {} # keep track of the largest pod modification recorded for a sequence

    # Dictionary for hairpin spliced sequences: hash(HairpinName)->sequence
    hairpinseq = load_splicedfasta(Path (args.mappingtables) / fn_splicedsequences);

    # loading isomiR type annotations
    isomir_annotations = {}
    if fn_otherannotations != None:
        isomir_annotations = load_otherannotations(Path (args.mappingtables) / fn_otherannotations);
        
    # loading isomiR snp annotations
    isomir_snps = {}
    if fn_snptable != None:
        isomir_snps = load_snps(Path (args.mappingtables) / fn_snptable);

    # load coordinates for meta coordinates
    metacoordinates = {}
    metacoordinates = load_metacoordinates (Path (args.mappingtables) / fn_metacoordinates);

    # Hash for isomiRs that are exclusive, ambiguous, or SNPs to miRNA-space:  hash(isomiRseq)->count.
    # Will only hold values that are expressed (for speed optimization).
    fastfrag_exclusive_counts = {}
    fastfrag_notexclusive_counts = {}
    fastfrag_snps_counts = {}
    fastfrag_all_counts = {} # counts of all uniq read sequences

    # Creates a two element list of [exclusive, notexclusive] dictionaries
    lookup_table_list = load_lookup_table(Path (args.mappingtables) / fn_lookuptable);

    # Hash for isomiRs that are exclusive to miRNA-space:  hash(isomiRseq)->count
    fastfrag_exclusive = lookup_table_list[0]

    # Hash for isomiRs that are not exclusive to miRNA-space:  hash(isomiRseq)->count
    fastfrag_notexclusive = lookup_table_list[1]

    # read in fastq file
    print('Reading in fastq file')
    stat_totalfragmentreads_exclusive = 0
    stat_totalfragmentreads_notexclusive = 0
    stat_totalfragmentreads_snps = 0
    trimmedfastqfile = open_read_gz (args.trimmedfastqfile, None);

    global stat_totalstartingreads
    line_number = 1
    for line in trimmedfastqfile:
        if line_number == 1 or line_number == 3:
            line_number += 1
            continue
        elif line_number == 4:
            line_number = 1
            continue
        else:
            # Legacy stuff used for reading all 4 lines of code in the fastq file
            # line_header = line.rstrip().decode("utf-8")
            # line_seq = trimmedfastqfile.rstrip().decode("utf-8")
            # line_misc = trimmedfastqfile.readline().rstrip().decode("utf-8")
            # line_phred = trimmedfastqfile.readline().rstrip().decode("utf-8")
            line_seq = line.rstrip().decode("utf-8")
            line_number += 1

        stat_totalstartingreads += 1

        if not VALID.match(line_seq):
            print("Error, exiting: Unrecognized character in fastq sequence " + line_seq + ", only [ATCGN]'s allowed",
                  file=sys.stderr)
            sys.exit(1)
        else:
           if line_seq in fastfrag_all_counts:
              fastfrag_all_counts[line_seq] += 1;
           else:
              fastfrag_all_counts[line_seq] = 1;
           
           
    print('Analyzing fastq file')
    for line_seq, seq_count in fastfrag_all_counts.items():
      potentialmod = line_seq;
      stretch3p = ''; # map stretches of A/T/C/G 3p post-modifications
      foundexact = False;
      
      while (len (potentialmod) >= min_fragment_len):
         is_exclusive = potentialmod in fastfrag_exclusive
         is_notexclusive = potentialmod in fastfrag_notexclusive
         if is_exclusive and is_notexclusive:
            print('Error, exiting: Fragment ' + potentialmod + ' should not exist in both the exclusive and ambiguous buckets', file=sys.stderr)
            sys.exit(1)
            
         if is_exclusive: 
            stat_totalfragmentreads_exclusive += seq_count;
            fastfrag_exclusive_counts[line_seq] = seq_count;
         elif is_notexclusive:
            stat_totalfragmentreads_notexclusive += seq_count;
            fastfrag_notexclusive_counts[line_seq] = seq_count;
            
         if is_exclusive or is_notexclusive: # if it's in the lookup table, regardless of exclusivity
            foundexact = True;
            fastfrag_hairpin_annotations[line_seq] = []
            fastfrag_largestPostMod[line_seq] = 0
               
            largestmod = add_annotations (fastfrag_hairpin_annotations[line_seq],
                             line_seq,
                             potentialmod,
                             hairpinseq)
           
            # keep track of largest post-modification found for a given sequence 
            if largestmod != 0:
               fastfrag_largestPostMod[line_seq] = largestmod;
               
            break; # once a best hit has been found, stop
        
         if stretch3p == '':
            stretch3p = potentialmod[-1:];
        
         if (not (potentialmod[-1:] == stretch3p)) or (stretch3p == 'N'):
            break;

         potentialmod = potentialmod[:-1] # remove last character
         
      if not foundexact: # lets look in SNP table if not found exactly
         if line_seq in isomir_snps:
            stat_totalfragmentreads_snps += seq_count;
            fastfrag_snps_counts[line_seq] = seq_count;
            
    trimmedfastqfile.close()

    # Create output files for both the exclusive and non-exclusive fragments
    print('Generating output:')
    exclusive_file = create_output(fastfrag_exclusive_counts,
                                   fastfrag_hairpin_annotations,
                                   isomir_annotations,
                                   fastfrag_largestPostMod,
                                   stat_totalfragmentreads_exclusive,
                                   'exclusive-' + DEFAULT_FRAGMENTTYPE + 's',
                                   metacoordinates,
                                   args.outputprefix,
                                   args.customrpm);
    ambiguous_file = create_output(fastfrag_notexclusive_counts,
                                   fastfrag_hairpin_annotations,
                                   isomir_annotations,
                                   fastfrag_largestPostMod,
                                   stat_totalfragmentreads_notexclusive,
                                   'ambiguous-' + DEFAULT_FRAGMENTTYPE + 's',
                                   metacoordinates,
                                   args.outputprefix,
                                   args.customrpm);
    if fn_snptable != None: # only create SNP output if the snp table is in the mappingbundle
       create_snp_output (fastfrag_snps_counts,
                          isomir_snps,
                          stat_totalfragmentreads_snps,
                         'snps-' + DEFAULT_FRAGMENTTYPE + 's',
                          args.outputprefix,
                          args.customrpm);

    print('Beginning conversions to GFF3:')

    create_gff3_output(args.outputprefix, [exclusive_file, ambiguous_file, '--meta', 'miRCarta', '--mt_name', mappingtables_name, '--mt_rev', mappingtables_rev, '--mt_description', mappingtables_description])
    create_gff3_output(args.outputprefix, [exclusive_file, ambiguous_file, '--meta', 'miRBase', '--mt_name', mappingtables_name, '--mt_rev', mappingtables_rev, '--mt_description', mappingtables_description])
    print('Finished')
