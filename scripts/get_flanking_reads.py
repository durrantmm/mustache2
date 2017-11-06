import sys
import pysam
import gzip

def main():
    genome_sam = sys.argv[1]
    insertseq_sam = sys.argv[2]
    fastq1 = sys.argv[3]
    fastq2 = sys.argv[4]
    class1 = sys.argv[5]
    class2 = sys.argv[6]
    genome_flanks_sam = sys.argv[7]
    genome_noflanks_sam = sys.argv[8]
    insertseq_flanks_sam = sys.argv[9]
    insertseq_noflanks_sam = sys.argv[10]

    get_flanking_reads(genome_sam, insertseq_sam, fastq1, fastq2, class1, class2, genome_flanks_sam,
                       genome_noflanks_sam, insertseq_flanks_sam, insertseq_noflanks_sam)


def read_class(path):
    with gzip.open(path) as infile:
        for line in infile:
            line = line.decode('UTF-8').strip().split()
            out = type("", (), {})()
            out.name = line[1]
            out.taxon = line[2]
            yield out


def read_sam_pairs(sam_path):

    while True:
        sam_file = pysam.AlignmentFile(sam_path, 'r')

        while True:
            try:
                out = type("", (), {})()
                out.r1 = next(sam_file)
                out.r2 = next(sam_file)
                yield out
            except StopIteration:
                break


def reverse_complement(seq):
    out = ''
    for c in reversed(seq):
        if c == 'A':
            out += 'T'
        elif c == 'T':
            out += 'A'
        elif c == 'G':
            out += 'C'
        elif c == 'C':
            out += 'G'
        else:
            out += c
    return out


def add_taxon(r1, r2, class1, class2):
    if r1.is_read1:
        r1.set_tag('TX', class1.taxon)
        r2.set_tag('TX', class2.taxon)
        out_gr = [r1, r2]
    else:
        r1.set_tag('TX', class2.taxon)
        r2.set_tag('TX', class1.taxon)
        out_gr = [r2, r1]

    return out_gr


def adjust_read_order(gr1, gr2, ir1, ir2, class1, class2):

    if gr1.is_read1:
        gr1.set_tag('TX', class1.taxon)
        gr2.set_tag('TX', class2.taxon)
        out_gr = [gr1, gr2]
    else:
        gr1.set_tag('TX', class2.taxon)
        gr2.set_tag('TX', class1.taxon)
        out_gr = [gr2, gr1]

    if ir1.is_read1:
        ir1.set_tag('TX', class1.taxon)
        ir2.set_tag('TX', class2.taxon)
        out_ir = [ir1, ir2]
    else:
        ir1.set_tag('TX', class2.taxon)
        ir2.set_tag('TX', class1.taxon)
        out_ir = [ir2, ir1]

    return out_gr+out_ir

def get_flanking_reads(genome_sam_path, insertseq_sam_path, fastq1_path, fastq2_path, class1_path, class2_path,
                       genome_flanks_sam_path, genome_noflanks_sam_path, insertseq_flanks_sam_path,
                       insertseq_noflanks_sam_path):

    genome_template = pysam.AlignmentFile(genome_sam_path, 'r')
    insertseq_template= pysam.AlignmentFile(insertseq_sam_path, 'r')

    genome_flanks_sam = pysam.AlignmentFile(genome_flanks_sam_path, 'w', template=genome_template)
    genome_noflanks_sam = pysam.AlignmentFile(genome_noflanks_sam_path, 'w', template=genome_template)
    insertseq_flanks_sam = pysam.AlignmentFile(insertseq_flanks_sam_path, 'w', template=insertseq_template)
    insertseq_noflanks_sam= pysam.AlignmentFile(insertseq_noflanks_sam_path, 'w', template=insertseq_template)

    genome_sam = read_sam_pairs(genome_sam_path)
    insertseq_sam = read_sam_pairs(insertseq_sam_path)

    class1 = read_class(class1_path)
    class2 = read_class(class2_path)

    genome_reads = next(genome_sam)
    insertseq_reads = next(insertseq_sam)
    count = 0

    print("Beginning Read Processing")
    for c1, c2 in zip(class1, class2):
        count += 1

        if c1.name == genome_reads.r1.query_name == insertseq_reads.r1.query_name:

            gr1, gr2, ir1, ir2 = adjust_read_order(genome_reads.r1, genome_reads.r2,
                                                   insertseq_reads.r1, insertseq_reads.r2,
                                                   c1, c2)
            gr1.set_tag('IS', ir2.reference_name)
            gr2.set_tag('IS', ir1.reference_name)
            ir1.set_tag('GN', gr2.reference_name)
            ir2.set_tag('GN', gr1.reference_name)

            gr1_mapped = not gr1.is_unmapped
            gr2_mapped = not gr2.is_unmapped
            ir1_mapped = not ir1.is_unmapped
            ir2_mapped = not ir2.is_unmapped

            if gr1_mapped and gr2_mapped and ir1_mapped and ir2_mapped:
                #print("WITHIN INSERTSEQ INSIDE GENOME")
                genome_noflanks_sam.write(gr1)
                genome_noflanks_sam.write(gr2)
                insertseq_noflanks_sam.write(ir1)
                insertseq_noflanks_sam.write(ir2)

            elif gr1_mapped and (not gr2_mapped) and (not ir1_mapped) and (ir2_mapped):
                #print("LEFT FLANK ABSENT INSERTSEQ")
                genome_flanks_sam.write(gr1)
                insertseq_flanks_sam.write(ir2)

            elif (not gr1_mapped) and gr2_mapped and ir1_mapped and (not ir2_mapped):
                #print("RIGHT FLANK ABSENT INSERTSEQ")
                genome_flanks_sam.write(gr2)
                insertseq_flanks_sam.write(ir1)

            elif gr1_mapped and gr2_mapped and (not ir1_mapped) and ir2_mapped:
                #print("LEFT FLANK PRESENT INSERTSEQ")
                genome_flanks_sam.write(gr1)
                insertseq_flanks_sam.write(ir2)

            elif gr1_mapped and gr2_mapped and ir1_mapped and (not ir2_mapped):
                genome_flanks_sam.write(gr2)
                insertseq_flanks_sam.write(ir1)

            genome_reads = next(genome_sam)
            insertseq_reads = next(insertseq_sam)

        elif c1.name == genome_reads.r1.query_name:
            gr1, gr2 = add_taxon(genome_reads.r1, genome_reads.r2, c1, c2)
            genome_noflanks_sam.write(gr1)
            genome_noflanks_sam.write(gr2)
            genome_reads = next(genome_sam)

        elif c1.name == insertseq_reads.r1.query_name:
            ir1, ir2 = add_taxon(insertseq_reads.r1, insertseq_reads.r2, c1, c2)
            insertseq_noflanks_sam.write(ir1)
            insertseq_noflanks_sam.write(ir2)
            insertseq_reads = next(insertseq_sam)

    print("Reads Processed:", count)

if __name__ == '__main__':
    main()