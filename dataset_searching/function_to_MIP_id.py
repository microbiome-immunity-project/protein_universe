#! /usr/bin/env python

import sys
import csv
import json
import argparse
import glob
import gzip
import multiprocessing
from collections import defaultdict

def get_mip_ids_function_scores_from_file( filename ) :
    '''Open file and returns lists of mip ids and function scores'''
    with gzip.open( filename, "rb" ) as json_handle:
        print( f"LOADING {filename}")
        json_data = json.load( json_handle )
        return json_data[ 'pdb_chains' ], json_data[ 'Y_hat' ]

def get_function_ids_function_desc_from_file( filename ) :
    '''Opens file and returns lists of function ids and function descriptions'''
    with gzip.open( filename, "rb" ) as json_handle:
        json_data = json.load( json_handle )
        return json_data[ 'goterms' ], json_data[ 'gonames' ]

#def search_for_functions_in_file( function_list, threshold, filename ) :
def search_for_functions_in_file( vt ) :
    '''Search for listed function in file, returns a dict of results that exceed threshold'''
    function_list, threshold, filename = vt
    #print(function_list, threshold, filename)
    mip_ids, scores = get_mip_ids_function_scores_from_file( filename )
    mip_ids_scores_found = defaultdict( list )

    for function in function_list:
        ontology, function_name, function_index = function
        for i, score_vec in enumerate( scores ):
            if scores[i][ function_index ] >= threshold:
                mip_ids_scores_found[ function ].append( ( mip_ids[i], scores[i][ function_index ] ) )

    return mip_ids_scores_found

def main():
    # arguments
    parser = argparse.ArgumentParser( description="A script for searching throught the DeepFRI outputs of MIP models. Loading all ontologies requires signifigant amounts of memory. The most efficient usage is to submit multipe functions at once." )

    parser.add_argument("-f", "--functions", type=str, nargs='+', help="List of GO or EC terms to look up seperated by whitespace. e.g. GO:0030246 EC:4.99.1.-", required=True )
    parser.add_argument("-t", "--threshold", type=float, help="Threshold for minimum DeepFRI score, range 0.0--1.0", default=0.10 )
    parser.add_argument("-d", "--deepfri_functions", type=str, help="Path to gziped json formatted DeepFRI output files from the supporting dataset, eg. rosetta_high_quality_function_predictions", required=True)
    parser.add_argument("-p", "--output_prefix", type=str, help="Prefix for name of output files", default="MIP_FUNCTIONS")
    parser.add_argument("-n", "--num_proc", type=int, help="Number of processes to run", default="8")


    args = parser.parse_args()

    # get function dicts
    ontologies = [ "BP", "MF", "CC", "EC" ]

    # dicts to hold everything
    function_ids = {}
    function_desc = {}

    # get function list ( the list of functions in each DeepFRI output is in the same order )
    for ontology in ontologies:
        function_ids[ ontology ], function_desc[ ontology ] = get_function_ids_function_desc_from_file( f"{ args.deepfri_functions }/DeepFRI_MIP_00000000_{ ontology }_pred_scores.json.gz" )

    # zip function, index, ontology for all ontology
    function_tups = []
    for ontology in ontologies:
        for i, id in enumerate( function_ids[ ontology ] ) :
            function_tups.append( ( ontology, id, i ) )

    # find ontology and index, for each function to search, strip 'EC:' from name
    function_search_data = []
    ontologies_to_search = set()
    for function in args.functions :
        if function[:3] == "EC:":
            function = function[3:]

        for ft in function_tups:
            if ft[1] == function:
                ontologies_to_search.add( ft[0] )
                function_search_data.append( ft )
                break


    # get list of all mf, bp, cc, ec files
    file_lists = {}
    for ontology in ontologies_to_search:
        file_lists[ ontology ] = glob.glob( f"{ args.deepfri_functions }/DeepFRI_MIP_*_{ ontology }_pred_scores.json.gz" )


    # search through function data and print out results
    for ontology in ontologies_to_search:

        # get lists of functions for this ontology and zip it up in to a single argument to pass to the multiprocess worker function
        functions_ontology = [ i for i in function_search_data if i[0] == ontology ]
        mp_args = zip( [ functions_ontology ]*len(file_lists[ ontology ] ), [args.threshold]*len(file_lists[ ontology ]), file_lists[ ontology ] )

        # search each file
        for func in functions_ontology:
            print( f"SEARCHING for { func[1] } in {func[0]} ontology" )

        pool = multiprocessing.Pool( processes=args.num_proc )
        function_data = pool.map( search_for_functions_in_file, mp_args )

        # output results for each function
        for func in functions_ontology:
            ont, f_idn, f_idx = func
            output_filename=f"{args.output_prefix}_{ont}_{f_idn}_{args.threshold}.csv"
            with open( output_filename, "w" ) as csvfile:
                csv_writer = csv.writer( csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                csv_writer.writerow( [ "MIP_ID", "DeepFRI_Score", "GO/EC_ID", "GO/EC_desc"])

                for file_result in function_data:
                    if func in file_result:
                        for result in file_result[func]:
                            csv_writer.writerow( [
                                result[0], # MIP id
                                result[1], # deepfri score
                                function_ids[ ont ][ f_idx ], # GO/EC id
                                function_desc[ ont ][ f_idx ] # GO/EC name
                            ] )


if __name__ == "__main__":
    main()
