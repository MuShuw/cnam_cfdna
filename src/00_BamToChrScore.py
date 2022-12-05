#!/usr/bin/env python3
import argparse
from multiprocessing import Pool
import numpy as np
import os
import pysam
import sys
import time


class ChromProcessor(object):
    """
    Class used for parallel processing of chromosomes.
    """

    def __init__(self, parameters):
        """
        Create the parameter for the chromosomes processing.

        Parameters
        ----------
        parameters : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.parameters = parameters

    def __call__(self, chromNsize):
        """
        Run the function for the chromosomes processing.

        Parameters
        ----------
        chromNsize : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        chr, size = chromNsize
        data_print(chr, size, self.parameters)
        return


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system and if
    it is a bam file
    -----
    Parameters :
        parser : argparse object
        arg : str
    -----
    Returns :
        file
    """
    file = os.path.abspath(arg)
    if not os.path.exists(file):
        parser.error("The file {} does not exist !".format(file))
    elif file.split("/")[-1].split(".")[-1] != "bam":
        parser.error("The file {} is not a bam file !".format(file))
    else:
        return file


def is_odd_number(w_size):
    """
    Check if the windows size is an odd number
    -----
    Parameters :
        w_size : (int) size of the windows
    -----
    Returns :
        w_size : (int) size of the windows
    """
    if w_size % 2 == 0:
        print("the window size has to be odd!")
        sys.exit()
    else:
        return w_size


def get_parser():
    """
    Get parser object for the script.
    -----
    Parameters :
        None
    -----
    Returns :
        parser : parser that return all the information written in the command
        Shell
    """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input",
                        dest="input",
                        type=lambda x: is_valid_file(parser, x),
                        help="input file to enter, the file has to be a bam file",
                        required=True)

    parser.add_argument("-s", "--wsize",
                        dest="wind_size",
                        default=121,
                        type=int,
                        help="size of the windows",
                        required=False)

    parser.add_argument("-l", "--frag_min",
                        dest="frag_min",
                        default=120,
                        type=int,
                        required=False,
                        help="minimum fragment size to keep")

    parser.add_argument("-m", "--frag_max",
                        dest="frag_max",
                        default=180,
                        type=int,
                        help="maximum fragment size to keep",
                        required=False)

    parser.add_argument("-d", "--sampling",
                        dest="sampling",
                        default=False,
                        help="activate the downsampling",
                        required=False)
    parser.add_argument("-", "--nb_core",
                        dest="nb_core",
                        default=10,
                        type=int,
                        help="choose the number of core to be used",
                        required=False)
    return parser


def get_chromosomes_names(input):
    """
    get the chromosome name and size of each chromosome using the header of the bam
    file
    -----
    Parameters :
        input : (str) name of the input file
    -----
    Returns :
        list_chromosomes : (list) list of all the chromosomes avaiblable present in the file
        list_length : (list) list of all the chromosomes size avaiblable  in the file
    """

    # opening the bam file with pysam
    bamfile = pysam.AlignmentFile(input, 'rb')
    # query all the names of  the chromosomes in a list
    list_chromosomes = bamfile.references[0:24]
    list_length = bamfile.lengths[0:24]
    bamfile.close()
    return list_chromosomes, list_length


def count_depth(chr, size, input):
    """
    Count the read depth.

    For each genomic coordinate return the number
    of reads the array return is a zero base array which means that the index
    begins with 0 instead of 1 for normal genomic coordinates
    -----
    Parameters :
        chr : (str) name of the chromosome
        size : (int) number of base in the chromosome
        input : (str) name of the input file
    -----
    Returns :
        array_depth : (numpy array) : depth of the file for each genomic
        coordinates
    """
    bamfile = pysam.AlignmentFile(input, 'rb')
    # WARNING this might reduce the range of the coverage
    array_depth = np.zeros(size+1, dtype=np.uint16)
    for pileupcolumn in bamfile.pileup(chr, stepper='nofilter'):
        pos = pileupcolumn.reference_pos
        depth = pileupcolumn.nsegments
        array_depth[pos+1] = depth
    bamfile.close()
    return array_depth


def get_wps_and_frag(chr, input, size, w_size, frag_min, frag_max):
    """
     Compute the WPS(Windows Protection Score) and  Get the size, start and
     stop coordinates for each fragment
    -----
    Parameters :
        chr : (str) name of the chromosome
        size : (int) number of base in the chromosome
        input : (str) name of the input file
        w_size : (int) size of the window
        frag_min: (int) minimal size for the fragment to be retain for wps
        frag_max: (int) maximal size for the fragment to be retain for wps
    -----
    Returns :
        array_overlap: (numpy array) : fragment depth or coverage of the file
        for each genomic coordinates
        array_wps: (numpy array) : fragment depth or coverage of the file
        for each genomic coordinates
    """
    bamfile = pysam.AlignmentFile(input, 'rb')
    # WARNING this might reduce the range of the coverage
    frag_overlap = np.zeros(size+1, dtype=np.uint16)
    wps_array = np.zeros(size+1, dtype=np.int16)
    for read in bamfile.fetch(contig=chr):
        if read.is_proper_pair:
            if read.mate_is_unmapped:
                continue
            if read.next_reference_id != read.reference_id:
                continue
            flen = read.template_length
            if frag_min <= flen <= frag_max:
                read_start = read.reference_start
                rlen = read.query_length
                read_end = read.reference_end
                next_read_start = read.next_reference_start
                next_read_end = next_read_start + rlen-1
                list_read = (read_start, read_end,
                             next_read_start, next_read_end)
                frag_start = min(list_read)
                frag_end = max(list_read)
                frag_overlap[frag_start:frag_end] += 1
                wps_array[int(frag_start - (w_size-1)/2):
                          int(frag_start+flen+(w_size-1)/2+1)] -= 1
                if (flen >= w_size):
                    wps_array[int(frag_start + (w_size-1)/2):
                              int(frag_start+flen-(w_size-1)/2)+1] += 2
    bamfile.close()
    return frag_overlap, wps_array


def data_print(chr, size, parameters):
    """
     Gathering of all the score compute in order to prepare an output with a
     TSV file format.

    -----
    Parameters :
        chr : (str) name of the chromosome
        size : (int) number of base in the chromosome
        input : (str) name of the input file
        w_size : (int) size of the window
        frag_min: (int) minimal size for the fragment to be retain
        frag_max: (int) maximal size for the fragment to be retain
    -----
    Returns :

    """
    input, w_size, frag_min, frag_max = parameters
    print(chr + "  is processing")
    out_name = input.split("/")[-1].split(".")[0] + "_" + chr

    depth_array = count_depth(chr, size, input)
    frag_array, wps_array = get_wps_and_frag(
        chr, input, size, w_size, frag_min, frag_max)
    t1 = time.time()
    depth_array, frag_array, wps_array = depth_array[1:, ],
    frag_array[1:, ], wps_array[1:, ]
    out_matrix = np.rec.fromarrays(
        [depth_array, frag_array, wps_array], names='rd,fd,wps')
    np.savez_compressed(out_name, out_matrix)
    t2 = time.time()
    total_time = round((t2 - t1)/60)
    print("total time for {} is {} min".format(chr, total_time))


def get_stat(input):
    """
     give for each chromosome the number of mapped and unmapped reads
    -----
    Parameters :
        input : (str) name of the input file
    -----
    Returns :
        stat: (list) : mapped and unmapped reads
    """
    bamfile = pysam.AlignmentFile(input, 'rb')
    stat = bamfile.get_index_statistics()
    bamfile.close()
    return stat


if __name__ == "__main__":
    args = get_parser().parse_args()
    print("input file : "+args.input)
    print("windows size is "+str(args.wind_size)+" bp")
    print("fragments used for WPS is between {} and {} bp".format(
        str(args.frag_min), str(args.frag_max)))
    print("number of core used : {}".format(str(args.nb_core)))
    list_chromosomes, list_length = get_chromosomes_names(args.input)
    is_odd_number(args.wind_size)
    # marche avec un parametres
    # test avec deux parametres
    try:
        #pool = Pool(os.cpu_count() - 1)
        pool = Pool(23)
        dataProcessor = ChromProcessor(
            (args.input, args.wind_size, args.frag_min, args.frag_max))
        pool.map(dataProcessor, zip(list_chromosomes, list_length))
    finally:
        pool.close()
        pool.join()
