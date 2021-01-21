To download the latest release of isoMiRmap please visit [here](https://cm.jefferson.edu/isoMiRmap/) or the 'Releases' section of this GitHub page.

# isoMiRmap
IsoMiRmap is an open source python application for the fast, deterministic, and exhaustive mining of isomiRs. The tool produces output in HTML, miRGFF3, and in tab-delimited formats.  

- Click [here](HOWTO_and_LICENSE.txt) to view the help and license for the tool
- Click [here](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btab016/35927328/btab016.pdf) to view the isoMiRmap publication.
- Click [here](https://cm.jefferson.edu/isoMiRmap/) to view more program details or to access downloadable isomiR profiles for 11,085 datasets from the 1000 Genomes and TCGA projects consisting of 5 population groups and 33 cancer types respectively.

# Abstract (from publication above)
### Motivation: 
MicroRNA (miRNA) precursor arms give rise to multiple isoforms simultaneously called “isomiRs.” IsomiRs from the
same arm typically differ by a few nucleotides at either their 5´ or 3´ termini, or both. In humans, the identities and abundances of
isomiRs depend on a person’s sex, population of origin, race/ethnicity, and on tissue type, tissue state, and disease type/subtype.
Moreover, nearly half of the time the most abundant isomiR differs from the miRNA sequence found in public databases. Accurate
mining of isomiRs from deep sequencing data is thus important.

### Results: 
We developed isoMiRmap, a fast, standalone, user-friendly mining tool that identifies and quantifies all isomiRs by directly
processing short RNA-seq datasets. IsoMiRmap is a portable “plug-and-play” tool, requires minimal setup, has modest computing
and storage requirements, and can process an RNA-seq dataset with 50 million reads in just a few minutes on an average laptop.
IsoMiRmap deterministically and exhaustively reports all isomiRs in a given deep sequencing dataset and quantifies them accurately
(no double-counting). IsoMiRmap comprehensively reports all miRNA precursor locations from which an isomiR may be transcribed,
tags as ‘ambiguous’ isomiRs whose sequences exist both inside and outside of the space of known miRNA sequences and reports the
public identifiers of common single-nucleotide polymorphisms and documented somatic mutations that may be present in an isomiR.
IsoMiRmap also identifies isomiRs with 3´ non-templated post-transcriptional additions. Compared to similar tools, isoMiRmap is
the fastest, reports more bona fide isomiRs, and provides the most comprehensive information related to an isomiR’s transcriptional
origin.

# Computational Medicine Center at Thomas Jefferson University
Please view our other research [here](https://cm.jefferson.edu/).  
