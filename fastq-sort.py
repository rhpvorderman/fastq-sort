import argparse
import functools
from operator import attrgetter
from typing import List
import heapq
import contextlib
import os

import dnaio
from xopen import xopen

READ_OPENER = functools.partial(xopen, threads=0)


def fastq_to_sorted_chunks(fastq: str, output_prefix: str,
                           max_ram_records: int = 5_000_000) -> List[str]:
    sorted_chunks = []
    with dnaio.open(fastq, mode="r", opener=READ_OPENER) as fastq_reader:
        while True:
            filename = output_prefix + str(len(sorted_chunks))
            records = [record for i, record in
                       zip(range(max_ram_records), fastq_reader)]
            if len(records) == 0:
                break
            records.sort(key=attrgetter('sequence'))
            with xopen(filename, mode="wb", format="gz", threads=0, compresslevel=1
                       ) as sorted_chunk:
                for record in records:
                    sorted_chunk.write(record.fastq_bytes())
            # If not explicitly deleted python keeps the reference alive
            # throughout the next loop, effectively doubling the memory usage.
            del(records)
            sorted_chunks.append(filename)
    return sorted_chunks


def sorted_chunks_to_sorted_fastq(sorted_chunks: List[str], output_file: str):
    with contextlib.ExitStack() as stack:
        files = [stack.enter_context(
                 dnaio.open(sorted_chunk, mode="r", opener=READ_OPENER))
                 for sorted_chunk in sorted_chunks]
        with xopen(output_file, mode="wb", compresslevel=1, threads=0
                   ) as output:
            for record in heapq.merge(*files, key=attrgetter('sequence')):
                output.write(record.fastq_bytes())


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", metavar="INPUT", help="Input FASTQ file")
    parser.add_argument("-o", "--output", help="Output FASTQ file")
    parser.add_argument("--tmp", help="Temporary sorted files are stored here")
    return parser


def main():
    args = argument_parser().parse_args()
    output_prefix = args.output
    if args.tmp:
        output_prefix = os.path.join(args.tmp, os.path.basename(args.output))
    sorted_chunks = fastq_to_sorted_chunks(args.input, output_prefix)
    sorted_chunks_to_sorted_fastq(sorted_chunks, args.output)


if __name__ == "__main__":
    main()

