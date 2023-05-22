#! /usr/bin/env python3
from Bio import SeqIO
from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call
from math import ceil
import rich_click as click
import subprocess
import os
import gzip
import random

@click.command()
@click.option('-a', '--assembly', type=str)
@click.option('-f', '--fastq', type=str)
@click.option('-m', '--model', type=str)
@click.option('-p', '--prefix', type=str)
@click.option('-t', '--threads', type=str)

def main(assembly, fastq, model, prefix, threads):
    if not os.path.exists(os.path.basedir(prefix)):
        os.mkdir(os.path.basedir(prefix))
    mini_align(fastq, assembly, prefix, threads)
    medaka_consensus(assembly, prefix, model, threads)
    medaka_stich(threads, assembly, prefix)

def mini_align(fastq, assembly, prefix, threads):
    if os.path.exists(assembly + ".mmi") and os.path.exists(prefix + ".bam"):
        print("Index and bam file are already present")
        return
    arguments = ["mini_align", "-i", fastq, "-r", assembly, "-m", "-p", prefix, "-t", threads]
    command = ' '.join(arguments)
    print(">", command)
    subprocess.run(command, shell=True, check=True)

def medaka_consensus(assembly, prefix, model, threads):
    # Read genome file to obtain contig information
    if (assembly.endswith(".gz")):
        reader = gzip.open(assembly, 'rb')
    else:
        reader = open(assembly)

    recids = []
    for rec in SeqIO.parse(reader, 'fasta'):
        recids.append(rec.id)
    batchnum = min(25, ceil(int(threads)/2))
    batchsize = int(len(recids)/batchnum + 1)
    reader.close()

    # Run this in parallel for each chunk by number of threads / 2
    commands = []
    for i in range(1,batchnum):
        batch = "{}_batch{}.hdf".format(prefix, i)
        contigs = random.choices(recids, k=batchsize)
        # Quiet due to the parallel nature and the many processes
        arguments = ["medaka", "consensus", "--quiet", prefix + ".bam", batch, "--model", model, "--threads", "1", "--region", ' '.join(contigs)]
        arguments = ' '.join(map(str,arguments))
        commands.append(arguments)
    
    # Multi thread section
    print("Number of jobs to be executed in parallel:", len(commands), "with maximum of", batchnum, "jobs")
    pool = Pool(batchnum) # two concurrent commands at a time
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("command failed with exit code %d and command %s" % (returncode, commands[i]))
        else:
            print("command finished: %d %s" % (returncode, commands[i]))

def medaka_stich(threads, assembly, prefix):
    arguments = ["medaka", "stitch", "--threads", threads, prefix + "*.hdf", assembly, "consensus.fasta"]
    command = ' '.join(arguments)
    print(">", command)
    subprocess.run(command, shell=True, check=True)


if __name__ == '__main__':
    main()
