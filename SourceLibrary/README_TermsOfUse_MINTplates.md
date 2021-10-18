# MINTplates "license-plates"

License Plates for tRNA-derived fragments (tRFs) and miRNA isoforms (isomiRs). Create sequence based names for your fragments.

This is the python License Plate (or MINTplate) code, please visit 
https://cm.jefferson.edu/MINTcodes/ to get the latest version.

## General Information

This code was created by Venetia Pliatsika, Isidore Rigoutsos, Jeffery Ma, Phillipe Loher
It can be used to create/encode molecular "license-plates" from sequences and to also decode the "license-plates"
back to sequences.  While initially created for tRFs (tRNA fragments), this tool can be used for 
any genomic sequences including but not limited to:  tRFs, isomiRs, reference miRNA, etc.
For more information on "license-plates", visit https://cm.jefferson.edu/MINTcodes/ and 
refer to publications https://www.ncbi.nlm.nih.gov/pubmed/27153631/ and https://www.ncbi.nlm.nih.gov/pubmed/28220888/.
Contact us at: https://cm.jefferson.edu/contact-us/

## License and Terms of Use

This MINTplates package is available under the open source GNU GPL v3.0 license
(https://www.gnu.org/licenses/gpl-3.0.en.html).

THE CODE IS PROVIDED “AS IS” WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESSED
OR IMPLIED. TO THE FULLEST EXTENT PERMISSIBLE PURSUANT TO APPLICABLE LAW. THOMAS JEFFERSON
UNIVERSITY, AND ITS AFFILIATES, DISCLAIM ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NON-INFRINGEMENT.

NEITHER THOMAS JEFFERSON UNIVERSITY NOR ITS AFFILIATES MAKE ANY REPRESENTATION AS TO THE RESULTS
TO BE OBTAINED FROM USE OF THE CODE.

## Usage Information

Usage (Python 2 and 3 compatible):

    python MINTplates.py example_sequences_to_encode.txt en --p [prefix to add to license plate]
    python MINTplates.py example_license_plates_to_decode.txt de --p [optional prefix, not used in decoding]
