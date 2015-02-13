#! /usr/bin/env python

"""
idability.py
============
Please type "./idability.py -h" for usage help

Authors:
  Eric A. Franzosa (franzosa@hsph.harvard.edu)
  Curtis Huttenhower (chuttenh@hsph.harvard.edu)

Copyright (c) 2015 Harvard T. H. Chan School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os, sys, argparse, csv

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = """
DESCRIPTION:

  This is a python program for generating and evaluating 
  hitting-set based codes. The program operates on tabular data 
  organized with samples  as columns and features as rows.
  (Also known as PCL format.)

BASIC OPERATION:

  When the program is given a table, it will attempt to construct 
  a unique set of features for each sample:

  $ ./idability.py demo1.pcl

  When given a table and a set of codes, the program will report 
  which samples are hit by which codes:

  $ ./idability.py demo2.pcl --codes demo1.codes.txt

  Without setting any additional arguments, the code construction
  process will be naive: only presence/absence information is
  considered and minimal codes are prioritized.

  Setting '--meta_mode [relab/rpkm]' will configure all settings
  to behave optimally for metagenomics features
  (measured in relative abundance [relab] or RPKM units)

  $ ./idability.py stool-markers-visit1.pcl --meta_mode rpkm

  is equivalent to:

  $ ./idability.py stool-markers-visit1.pcl -s 0.8 -m 7 -d 5 -n 0.05 -r abundance_gap

  Parameters can be fine-tuned for user-specific applications.

ARGUMENTS:
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_na = "#N/A"
c_epsilon = 1e-20
c_codes_extension = "codes.txt"
c_hits_extension = "hits.txt"

# ---------------------------------------------------------------
# arguments
# ---------------------------------------------------------------

def funcGetArgs ():
    """ master argument parser """

    parser = argparse.ArgumentParser( 
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument( 
        'table',
        type=str,
        help="""
        Tab-delimited table file to encode/decode.
        PCL format: rows are features, cols are samples (or subjects), both have headers
        """,
    )
    parser.add_argument( 
        "-c", '--codes',
        type=str,
        help="""
        Codes file produced in an earlier run.
        Specifying this option will compare codes to input table.
        """,
    )
    parser.add_argument( 
        "-j", '--jaccard_similarity_cutoff',  
        type=float,
        help="""
        If set, an encoded feature will 'knock out' similar features
        at the specified threshold (Jaccard score in [0-1]).
        """,
    )
    parser.add_argument( 
        "-m", '--min_code_size',
        type=int,
        default=1,
        help="""
        If set, codes will continue to be lengthened beyond
        the point of uniqueness. Limits spurious hits in time-varying data.
        """,
    )
    parser.add_argument( 
        "-d", '--abund_detect',
        type=float,
        default=c_epsilon,
        help="""
        Features with values above this are scored as confidently present.
        When running in encode mode, this restricts the features that can be added to a code.
        When running in decode mode, this restricts the features that can be hit by a code.
        """,
    )
    parser.add_argument( 
        "-n", '--abund_nondetect',        
        type=float,
        default=c_epsilon,
        help="""
        Features with values below this are scored as confidently absent.
        Only applied to encode mode. A sample with a feature below this threshold
        is considered to be 'missing' the feature for hitting set purposes.
        """,
    )
    parser.add_argument( 
        "-r", "--ranking",
        type=str,
        default="rarity",
        choices=["rarity", "abundance_gap"],
        help="""
        The method by which an individual's features should be prioritized when
        building codes. The default, rarity, prioritizess-less prevalent features.
        The alternative method, abundance_gap, prioritizes features with a large 
        abundance gap between the individual's value and the next highest value.
        """,
    )
    parser.add_argument( 
        "-o", "--output",
        type=str,
        help="""
        Name for the output file (codes or confusion matrix, depending on mode).
        If not supplied, a default will be constructed from the input file names.
        """,
    )
    parser.add_argument( 
        "-e", "--meta_mode",
        type=str,
        choices=["off", "relab", "rpkm"],
        default="off",
        help="""
        Automatically optimize all variables for working with metagenomic codes.
        If working with relative abundance data, select the "relab" mode.
        If working with reads per kilobase per million reads (RPKM) data, select the "rpkm" mode.
        """,
    )

    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# utilities and file i/o
# ---------------------------------------------------------------

def funcPathToName ( path ):
    return os.path.split( path )[1].split( "." )[0]

def funcLoadSFV ( path, cutoff ):
    """ 
    loads a table file to a nested dict (sfv: sample->feature->value)
    values below cutoff are ignored to save space
    """
    sfv = {}
    with open( path ) as fh:
        headers = None
        for row in csv.reader( fh, dialect="excel-tab" ):
            if headers is None:
                headers = row[1:]
                sfv = {header:{} for header in headers}
            else:
                feature, values = row[0], list(map( float, row[1:] ))
                assert len( values ) == len( headers ), \
                    "row length mismatch"
                for header, value in zip( headers, values ):
                    if value >= cutoff:
                        sfv[header][feature] = value
    return sfv

def funcReduceSFV ( sfv, cutoff, greater=True ):
    """
    rebuild sfv with only entries > cutoff
    maintain all samples(s), even if no features in sample meet cutoff
    """
    temp = {sample:{} for sample in sfv}
    for sample, fdict in list(sfv.items()):
        for feature, value in list(fdict.items()):
            if ( greater and value >= cutoff ) or ( not greater and value < cutoff ):
                temp[sample][feature] = value
    return temp

def funcFlipSFV ( sfv ):
    """
    make a fsv object, i.e. feature->sample->value map
    """
    fsv = {}
    for sample, fdict in list(sfv.items()):
        for feature, value in list(fdict.items()):
            fsv.setdefault( feature, {} )[sample] = value
    return fsv

def funcSetForm ( nested_dict ):
    """ 
    reduces inner dict to key set when we're done with values
    """
    return {k:set( list(v) ) for k, v in list(nested_dict.items())}

def funcCheckHits ( sample_hits ):
    """
    produces confusion results by comparing keys to lists of hit samples
    """
    counts = {k:0 for k in "1|TP 3|FN+FP 2|TP+FP 4|FN 5|NA".split()}
    for sample, hits in list(sample_hits.items()):
        if hits is None:
            counts["5|NA"] += 1
        else:
            tp_hit = True if sample in hits else False
            fp_hit = True if len([sample2 for sample2 in hits if sample2 != sample]) > 0 else False
            if tp_hit:
                key = "2|TP+FP" if fp_hit else "1|TP"
            else:
                key = "3|FN+FP" if fp_hit else "4|FN"
            counts[key] += 1
    return counts

def funcWriteCodes ( sample_codes, path ):
    """
    write code sets to a text file
    """
    with open( path, "w" ) as fh:
        fh.write("#SAMPLE\tCODE\n")
        for sample in sorted( list(sample_codes) ):
            code = sample_codes[sample]
            items = [sample] 
            items += [c_na] if code is None else code
            fh.write("\t".join( items )+"\n")
    sys.stderr.write("wrote codes to:" + path + "\n")

def funcReadCodes ( path ):
    """
    read back in the codes written by funcWriteCodes
    """
    sample_codes = {}
    with open( path ) as fh:
        fh.readline() # headers
        for line in fh:
            items = line.strip().split( "\t" )
            sample, code = items[0], items[1:]
            sample_codes[sample] = code if c_na not in code else None
    return sample_codes

def funcWriteHits ( sample_hits, path ):
    """
    write hit results and summary to a text file
    """
    # compute confusion line
    confusion = funcCheckHits( sample_hits )
    with open( path, "w" ) as fh:
        for confusion_class in sorted( confusion ):
            count = confusion[confusion_class]
            fh.write("# %s: %d" % ( confusion_class, count ) + "\n")
        for sample in sorted( sample_hits ):
            hits = sample_hits[sample]
            items = [sample]
            if hits is None:
                items += ["no_code", c_na]
            else:
                items += ["matches" if len( hits ) > 0 else "no_matches"]
                items += hits
            fh.write("\t".join( items ) + "\n")
    sys.stderr.write("wrote hits to:" + path + "\n")

# ---------------------------------------------------------------------------
# encode part
# ---------------------------------------------------------------------------

def funcJaccard ( set1, set2 ):
    """ 
    jaccard similarity for two sets
    """
    count_union = len( set1.__or__( set2 ) )
    count_intersection = len( set1.__and__( set2 ) )
    return count_intersection / float( count_union )

def funcRankAbundGap( sfv, fsv, abund_nondetect ):
    """ 
    abundance gap sorting sfv features
    """
    sorted_features = {}
    for sample, fdict in list(sfv.items()):
        gaps = {}
        for feature, focal_value in list(fdict.items()):
            lesser_values = [abund_nondetect]
            lesser_values += [v for k, v in list(fsv[feature].items()) \
                              if v <= focal_value and k != sample]
            gaps[feature] = focal_value - max( lesser_values )
        sorted_features[sample] = sorted( 
            list(gaps), key=lambda feature: gaps[feature], )
    return sorted_features

def funcRankRarity( sfv, fsv, abund_nondetect ):
    """
    rarity sorting of sfv features
    """
    sorted_features = {}
    for sample, fdict in list(sfv.items()):
        sorted_features[sample] = sorted( 
            list(fdict), key=lambda feature: len( fsv[feature] ), 
            reverse=True, )
    return sorted_features

def funcMakeOneCode ( sample, ranked_features, sfv_sets, fsv_sets, \
                      similarity_cutoff, min_code_size ):
    """ 
    execute the idabilty algorithm for one sample 
    """
    features = ranked_features[:]
    other_samples = {sample2 for sample2 in sfv_sets if sample2 != sample}
    code = []
    while len( features ) > 0 and \
          ( len( other_samples ) > 0 or len( code ) < min_code_size ):
        feature = features.pop()
        code.append( feature )
        # restrict other samples
        old_count = len( other_samples )
        other_samples = other_samples.__and__( fsv_sets[feature] )
        new_count = len( other_samples )
        # forget current feature if it doesn't improve things
        # *** unless we've already knocked everyone out and are just lengthening code ***
        if old_count == new_count and old_count != 0:
            code.pop()
        # restrict remaining features to avoid similarity to best feature
        if similarity_cutoff is not None:
            features = list(filter( lambda feature2: \
                               funcJaccard( fsv_sets[feature], fsv_sets[feature2] ) < \
                               similarity_cutoff, features ))
    return code if len( other_samples ) == 0 else None

def funcEncode ( sfv, abund_detect, abund_nondetect, similarity_cutoff, min_code_size, ranking="rarity" ):
    """    
    run idability algorithm on all samples 
    """
    # flip sfv to fsv
    fsv = funcFlipSFV( sfv )
    # rebuild sfv with only features above abund threshold
    sfv = funcReduceSFV( sfv, cutoff=abund_detect )
    # prioritize features
    sys.stderr.write("performing requested feature ranking:" + ranking + "\n")
    funcRank = {"rarity":funcRankRarity, "abundance_gap":funcRankAbundGap}[ranking]
    sorted_features = funcRank( sfv, fsv, abund_nondetect )
    # simplify sfv and fsv to sets
    sfv_sets = funcSetForm( sfv )
    fsv_sets = funcSetForm( fsv )
    # make codes for each sample
    sample_codes = {}
    for i, sample in enumerate( list(sfv_sets) ):
        sample_codes[sample] = funcMakeOneCode( 
            sample, 
            sorted_features[sample],
            sfv_sets,
            fsv_sets,
            similarity_cutoff, 
            min_code_size,
        )
    return sample_codes

# ---------------------------------------------------------------------------
# decode part
# ---------------------------------------------------------------------------

def funcCheckOneCode ( code, sfv_sets ):
    """ 
    compare a single code to a population
    """
    code_set = set( code )
    hits = []
    for sample, features_set in list(sfv_sets.items()):
        if code_set.issubset( features_set ):
            hits.append( sample )
    return hits

def funcDecode ( sfv, sample_codes, abund_detect ):
    """
    compare all codes to a population
    """
    sfv_sets = funcSetForm( funcReduceSFV( sfv, abund_detect ) )
    sample_hits = {}
    for sample, code in list(sample_codes.items()):
        sample_hits[sample] = None if code is None else funcCheckOneCode( code, sfv_sets )
    return sample_hits

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    
    """ """

    # process arguments
    args = funcGetArgs()
    table_path = args.table
    codes_path = args.codes
    abund_detect = args.abund_detect
    abund_nondetect = args.abund_nondetect
    similarity_cutoff = args.jaccard_similarity_cutoff
    min_code_size = args.min_code_size
    ranking = args.ranking
    output_path = args.output

    # overrides
    if args.meta_mode != "off":
        choice = args.meta_mode
        abund_detect = 5.0 if choice == "rpkm" else 0.001
        abund_nondetect = abund_detect / 100.0
        # relax detection parameter in decoding step
        abund_detect = abund_detect / 10.0 if args.codes is not None else abund_detect
        similarity_cutoff = 0.8
        min_code_size = 7
        ranking = "abundance_gap"

    # determine output file name
    if output_path is None:
        items = [funcPathToName( table_path )]
        if codes_path is None:
            items.append( c_codes_extension )
        else:
            items.append( funcPathToName( codes_path ) )
            items.append( c_hits_extension )
        output_path = ".".join( items )

    # do this for either encoding/decoding
    sys.stderr.write("loading table file:" + table_path + "\n")
    sfv = funcLoadSFV( table_path, abund_nondetect )

    # make codes mode
    if codes_path is None:
        sys.stderr.write("encoding the table\n")
        sample_codes = funcEncode( 
            sfv, 
            abund_detect=abund_detect, 
            abund_nondetect=abund_nondetect, 
            similarity_cutoff=similarity_cutoff,
            min_code_size=min_code_size,
            ranking=ranking,
        )
        funcWriteCodes( sample_codes, output_path )

    # compare codes to table mode
    else:
        sys.stderr.write("decoding the table\n")
        sample_codes = funcReadCodes( codes_path )
        sample_hits = funcDecode( 
            sfv, 
            sample_codes,
            abund_detect=abund_detect,
        )
        funcWriteHits( sample_hits, output_path )

if __name__ == "__main__":
    main()
