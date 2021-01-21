IsoMiRmap (available at https://cm.jefferson.edu/isoMiRmap/)
------------------------------

1. General Information
----------------------
This code was created by Phillipe Loher, Isidore Rigoutsos, and Jeffery Ma.  
IsoMiRmap can be used to generate isomiR profiles from Short-RNA datasets.

More information can be found at https://cm.jefferson.edu/isoMiRmap/

If used, please cite:
TBA


2. License
---------------
IsoMiRmap is available under the open source GNU GPL v3.0 license (https://www.gnu.org/licenses/gpl-3.0.en.html).  
MINTplates is included and is governed by a seperate license.  See SourceLibrary/README_TermsOfUse_MINTplates.txt for more information.


3. Other software requirements
---------------
Python 3.5 (or greater) is required.


4. Usage Information
--------------------

usage: IsoMiRmap.py [-h] [--m MAPPINGTABLES] [--p OUTPUTPREFIX]
                    [--d CUSTOMRPM]
                    trimmedfastqfile

Use this program to generate isomiR profiles. Requires Python 3.5 or later.
OUTPUT:
- Plain text and HTML file pairs for exclusive isomiRs (*.exclusive-isomiRs.expression.*). IsomiRs that exist exclusively within miR-space (higher c
onfidence) will be reported in this file.  Also, their counterparts with 3p additions (e.g. uridylation) will be included.  RPM and annotation information included.
- Plain text and HTML file pairs for ambiguous isomiRs (*.ambiguous-isomiRs.expression.*). Same as above but for isomiRs that exist both within miR-space and also elsewhere in the genome.
- Variant-containing isomiRs (*.snps-isomiRs.expression.*)
- High level mapping stats are also generated seperately for exclusive, non-exclusive, and variant containing isomiRs (*.countsmeta.txt)
- Output in mirGFF3 format (*.gff3)

required arguments:
  trimmedfastqfile      This is a required argument. The file
                        “trimmedfastqfile” contains the sequencer reads from
                        the Short-RNA dataset. Any trimming (e.g. quality and
                        adapter trimming) must be done prior to running this
                        tool. If “trimmedfastqfile” ends in “.gz” it will be
                        treated as a gzipped FASTQ file. The FASTQ file
                        contains four lines per read (more info here:
                        https://en.wikipedia.org/wiki/FASTQ_format). Color-
                        space reads are not supported.

optional arguments:
  -h, --help            show this help message and exit
  --m MAPPINGTABLES, --mappingbundle MAPPINGTABLES
                        If not specified, “mappingbundle” will be set to
                        "MappingBundles/miRCarta". This is the relative or
                        absolute path to the mapping bundle that will be used.
                        Out of the box, we provide mapping bundles for both
                        miRCarta and miRBase for Homo Sapiens.
  --p OUTPUTPREFIX, --outputprefix OUTPUTPREFIX
                        If not specified, “outputprefix” will be set to
                        "output" and all output files will be generated in the
                        current working directory. The “outputprefix” is a
                        string that will be prepended in the output files that
                        isoMiRmap_v5 will generate. The string is meant to
                        serve as a mnemonic for the user. “outputprefix” can
                        optionally include relative or absolute directory
                        paths (if the directories exist) to save the output to
                        a different directory.
  --d CUSTOMRPM, --customRPM CUSTOMRPM
                        This is an optional argument. The value “customRPM” is
                        meant to be used as alternative denominator when
                        computing the RPM abundance of an isomiR. When this
                        parameter is defined by the user, an additional column
                        in the output files will be populated with RPM values
                        that have been calculated using this value in the
                        denominator – i.e. these values will be equal to raw
                        reads/<customRPM>*1,000,000. A common value to use
                        here is the original number of sequenced reads prior
                        to quality- and adaptor-trimming.


5.  Example Run
--------------------
python IsoMiRmap.py Example/random10k.LCL_YRI_NA19257.fastq.gz
