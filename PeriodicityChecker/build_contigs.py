import subprocess


def build_contigs(fastq: str, tmp_dir: str, kmer_length: int, min_contig_length=100, 
                     num_threads=1, inchworm_path='inchworm') -> str:
    '''
    wrapper to run inchworm assembler on reads

    inputs:
        fastq: path to fastq file
        tmp_dir: path to tmp directory
        kmer_length: length of kmers to use
        min_contig_length: minimum contig length to report
        max_contig_length: maximum contig length to report
        num_threads: number of threads to use
        inchworm_path: path to inchworm executable  

    outputs:
        contigs: path to contigs fasta file
    '''
    contigs = tmp_dir + '/contigs.fa'
    cmd = [inchworm_path, '--reads', fastq, '-K', str(kmer_length),
              '-L', str(min_contig_length), '--num_threads', num_threads, '> ', contigs]

    subprocess.call(cmd)

    return contigs

