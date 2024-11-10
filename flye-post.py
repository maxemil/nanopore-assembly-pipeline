#! /usr/bin/env python3

from Bio import SeqIO,SeqUtils
import pandas as pd
import click

@click.command()
@click.argument('flye_dir', type=str)
@click.option('-m', '--min_length', type=float, default=1000)
@click.option('-p', '--prefix', type=str)

def main(flye_dir, min_length, prefix):
    seqdict = parse_sort_sequences(flye_dir)
    seqdict, old2new = rename_contigs(seqdict, prefix)
    df = rename_table(flye_dir, old2new, seqdict)
    df.to_csv('{}/{}.tsv'.format(flye_dir, prefix), sep='\t', header=True, index=False)
    write_sequences("{}/{}.fasta".format(flye_dir, prefix), seqdict, min_length)

def rename_contigs(seqdict, prefix):
    old2new = {}
    new_seqdict = {}
    for i, rec in enumerate(seqdict.values()):
        oldid = rec.id
        rec.id = "{}_{:0{digits}d}".format(prefix, i+1, digits=len(str(len(seqdict))))
        rec.description = "orig_id:{}".format(oldid)
        old2new[oldid] = rec.id
        new_seqdict[rec.id] = rec
    return new_seqdict, old2new

def rename_table(flye_dir, old2new, seqdict):
    df = pd.read_csv('{}/assembly_info.txt'.format(flye_dir), sep='\t')
    df['orig_id'] = df['#seq_name']
    df['#seq_name'] = df['#seq_name'].apply(lambda x: old2new[x])
    df['gc'] = df['#seq_name'].apply(lambda x: SeqUtils.gc_fraction(seqdict[x].seq))
    return df

def parse_sort_sequences(flye_dir):
    seqdict = {}
    for rec in SeqIO.parse('{}/assembly.fasta'.format(flye_dir), 'fasta'):
        seqdict[rec.id] = rec
    return {k: v for k, v in sorted(seqdict.items(), key=lambda item: len(item[1].seq), reverse=True)}

def write_sequences(outfile, seqdict, min_length ):
    with open(outfile, 'w') as out:
        for rec in seqdict.values():
            if len(rec.seq) >= min_length:
                SeqIO.write(rec, out, 'fasta')


if __name__ == '__main__':
    main()
