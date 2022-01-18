import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
import argparse


def fix_fasta(input, output, prefix=None):

    encoding = guess_type(input)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

    fixed_seqs = []
    with _open(input) as f:
        for record in SeqIO.parse(f, "fasta"):
            record.description = ""
            if bool(prefix):
                record.id = f"{prefix}_{record.id}"
            fixed_seqs.append(record)

    SeqIO.write(fixed_seqs, output, "fasta")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "Updates FASTA record ids and removes their descriptions."
    )
    parser.add_argument(
        "--input", metavar="FILE", help="input file in FASTA format", required=True
    )
    parser.add_argument("--output", metavar="FILE", help="output file", required=True)
    parser.add_argument(
        "--prefix", metavar="string prefix", help="String prepended to FASTA id"
    )
    args = parser.parse_args()

    fix_fasta(**args.__dict__)
