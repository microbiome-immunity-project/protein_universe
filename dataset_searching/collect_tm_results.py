#! /usr/bin/env python

'''
Parse a one to all TM-align log file and merge results them into a single CSV file with headers

chain_1_name, chain1_chain, chain1_model, chain1_length, chain_2_name, chain2_chain, chain2_model, chain2_length, aligned_length, aligned_length_rmsd, aligned_length_seq_ident, tm_score_norm_chain1, tm_score_norm_chain2
'''

import pickle
import argparse
import urllib.request
import csv
import glob
from collections import defaultdict

def parse_tm_output( tm_out_filename ):
    '''
    Parse the TM-align output file that looks like below. Return a list of dictonaries with keys as above

00  *********************************************************************
01  * TM-align (Version 20210224): protein structure alignment          *
02  * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
03  * Please email comments and suggestions to yangzhanglab@umich.edu   *
04  *********************************************************************
05
06 Name of Chain_1: ../ros_models_grad_ordered/MIP1_00245968_LE.pdb:1:A (to be superimposed onto Chain_2)
07 Name of Chain_2: 6J34.pdb:1:A
08 Length of Chain_1: 151 residues
09 Length of Chain_2: 930 residues
10
11 Aligned length= 92, RMSD=   4.77, Seq_ID=n_identical/n_aligned= 0.065
12 TM-score= 0.36347 (if normalized by length of Chain_1, i.e., LN=151, d0=4.58)
13 TM-score= 0.08379 (if normalized by length of Chain_2, i.e., LN=930, d0=10.24)
14 (You should use TM-score normalized by length of the reference structure)
15
16 (":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)
17 FDFNKNSDLSNW------TIVN-DVIMGGVSNSTLEINPDGNAVFSGTVSLENNGGFSLVRYRFETLKVQGYTKAIIRLKGDPTTYQFRAKTNSND-R---------YSYIAIFETTGE-------------WQTIEISLSEMYPAFRGRKLDLPNYQTEQMQ--EVAFLIGNKKEEKFELVI-------
18                   ..   .  : :..::::          :            ::..:.. . . .:::.::::::.::::::::: .:.. .         .::::::...               .:::::::::.             : :::::  :::::::::..:::::::
19 ------------DVVVRLVY--DS--R-ADAFRAA----------F------------GVALADA-H-W-VDKTTLLWPGGENKPIVRLYY-SHSSKVAADSNGEFSDKYVKLTPTT--VSQQVSMRFPHLASYPAFKLPDDV-------------N-VDELLQGETVAIAAESDGILSSATQVQTAGVL
    '''

    data = []

    with open( tm_out_filename, "r" ) as tm_out:
        lines = [ line.rstrip() for line in tm_out ]


    # were gonna jump around a bit so get an index
    i = 0
    while i < len(lines):

        # read until we get to a new entry denoted by " *****"
        if lines[ i ].startswith( " *****" ):

            # get chain{1,2}_{name,model,chain} from lines i+6 & i+7
            chain1_name, chain1_model, chain1_chain = lines[ i+6 ].split()[ 3 ].split( ":" )
            chain2_name, chain2_model, chain2_chain = lines[ i+7 ].split()[ 3 ].split( ":" )

            # get aligned_length, rmsd_aligned_length, sequence_id_aligned_length from lines i+11
            align_split = lines[ i+11 ].split()
            aligned_length = align_split[ 2 ].split( "," )[ 0 ]
            aligned_length_rmsd = align_split[ 4 ].split( "," )[ 0 ]
            aligned_length_seq_ident = align_split[ 6 ]

            # get tm_scor_norm_chain{1,2} from lines i+12 and i+13
            tm_score_norm_chain1 = lines[ i+12 ].split()[1]
            tm_score_norm_chain2 = lines[ i+13 ].split()[1]

            # add dict to list
            data.append( {
                "chain1_name" : chain1_name,
                "chain1_model" : chain1_model,
                "chain1_chain" : chain1_chain,
                "chain2_name" : chain2_name,
                "chain2_model" : chain2_model,
                "chain2_chain" : chain2_chain,
                "tm_score_norm_chain1" : tm_score_norm_chain1,
                "tm_score_norm_chain2" : tm_score_norm_chain2,
                "aligned_length" : aligned_length,
                "aligned_length_rmsd" : aligned_length_rmsd,
                "aligned_length_seq_ident" : aligned_length_seq_ident,
                "aligned_length" : aligned_length,
                "aligned_length_rmsd" : aligned_length_rmsd,
                "aligned_length_seq_ident" : aligned_length_seq_ident,
                } )

            # increment to next entry
            i += 19

        else:
            # increment to next line
            i += 1

    return data

if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tm_align_output", type=str, help="Path to tm-align output file")
    parser.add_argument("-c", "--merged_csv", type=str, help="Path to csv file of merged tm-align output files")

    args = parser.parse_args()

    # parse file, sort based on "tm_score_norm_chain1"
    tm_data = sorted( parse_tm_output( args.tm_align_output ), key=lambda d: d[ "tm_score_norm_chain1" ] )

    # write it out to a CSV file
    with open( args.merged_csv, "w" ) as csv_file:
        dict_writer = csv.DictWriter( csv_file, tm_data[0].keys() )
        dict_writer.writeheader()
        dict_writer.writerows( tm_data )
