from ._dehasher_impl import Dehasher

import argparse
from sourmash.plugins import CommandLinePlugin
from tqdm import tqdm
    
class Command_Dehasher(CommandLinePlugin):
    
    command = "dehash"
    description = "Dehash sourmash signatures to ACGT k-mers"
    formatter_class = argparse.RawTextHelpFormatter

    
    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument("--sig-paths", nargs="+", required=True, help="signature files to process")
        subparser.add_argument("--fasta-paths", nargs="+", required=True, help="FASTA files to process")
        subparser.add_argument('-k', "--ksize", type=int, required=True, help="k-mer size")
        subparser.add_argument('-c', "--cores", type=int, required=False, default = 1, help="Number of cores")
        subparser.add_argument('-o', "--out", type=str, required=True, help="output signature path")

    
    def main(self, args):
        
        
        print("Loading sourmash signatures ...")
        
        fasta_files = args.fasta_paths
        output = args.out
        
        dehasher_obj = Dehasher(
            kSize=args.ksize,
            sig_paths = args.sig_paths,
        )

        for fasta in fasta_files:
            print(f"Processing {fasta}...")
            if args.cores > 1:
                dehasher_obj.map_kmer_to_hashes_single_fasta_parallel(fasta, args.cores)
            else:
                dehasher_obj.map_kmer_to_hashes_single_fasta(fasta)
        
        dehasher_obj.dump_kmers_to_file(output)