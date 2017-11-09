import sys
import pysam
import gzip

def main():
    genome_sam = sys.argv[1]
    insertseq_sam = sys.argv[2]
    class1 = sys.argv[3]
    class2 = sys.argv[4]
    genome_flanks_sam = sys.argv[5]
    genome_flanks_mate_sam = sys.argv[6]
    genome_noflanks_sam = sys.argv[7]
    insertseq_flanks_sam = sys.argv[8]
    insertseq_flanks_mate_sam = sys.argv[9]
    insertseq_noflanks_sam = sys.argv[10]
    sample = sys.argv[11]

    get_flanking_reads(genome_sam, insertseq_sam, class1, class2, genome_flanks_sam, genome_flanks_mate_sam,
                       genome_noflanks_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam, insertseq_noflanks_sam, sample)


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

def set_flank_tags(gr1, gr2, ir1, ir2, filter, sample):
    
    gr1.set_tag('FL', filter)
    gr1.set_tag('IS', ir2.reference_name)
    gr1.set_tag('SP', sample)

    gr2.set_tag('FL', filter)
    gr2.set_tag('IS', ir1.reference_name)
    gr2.set_tag('SP', sample)

    ir1.set_tag('FL', filter)
    ir1.set_tag('GN', gr2.reference_name)
    ir1.set_tag('SP', sample)

    ir2.set_tag('FL', filter)
    ir2.set_tag('GN', gr1.reference_name)
    ir2.set_tag('SP', sample)
    
    return gr1, gr2, ir1, ir2

def write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam, ir1, insertseq_flanks_sam,
                            ir2, insertseq_flanks_mate_sam):

    genome_flanks_sam.write(gr1)
    genome_flanks_mate_sam.write(gr2)
    insertseq_flanks_sam.write(ir1)
    insertseq_flanks_mate_sam.write(ir2)

def is_soft_clipped_once(read):
    if (read.cigarstring[1] == 'S' and read.cigarstring[-1] != 'S'):
        return 1
    elif (read.cigarstring[1] != 'S' and read.cigarstring[-1] == 'S'):
        return 2
    return False

def is_soft_clipped_right(read):
    if read.cigarstring.count('S') == 1:
        return True
    return False

def get_flanking_reads(genome_sam_path, insertseq_sam_path, class1_path, class2_path,
                       genome_flanks_sam_path, genome_flanks_mate_sam_path, genome_noflanks_sam_path,
                       insertseq_flanks_sam_path, insertseq_flanks_mate_sam_path, insertseq_noflanks_sam_path,
                       sample="None"):

    genome_template = pysam.AlignmentFile(genome_sam_path, 'r')
    insertseq_template= pysam.AlignmentFile(insertseq_sam_path, 'r')

    genome_flanks_sam = pysam.AlignmentFile(genome_flanks_sam_path, 'w', template=genome_template)
    genome_flanks_mate_sam = pysam.AlignmentFile(genome_flanks_mate_sam_path, 'w', template=genome_template)
    genome_noflanks_sam = pysam.AlignmentFile(genome_noflanks_sam_path, 'w', template=genome_template)
    insertseq_flanks_sam = pysam.AlignmentFile(insertseq_flanks_sam_path, 'w', template=insertseq_template)
    insertseq_flanks_mate_sam = pysam.AlignmentFile(insertseq_flanks_mate_sam_path, 'w', template=insertseq_template)
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
            #print("GENOME-IS ALIGNMENT")

            gr1, gr2, ir1, ir2 = adjust_read_order(genome_reads.r1, genome_reads.r2,
                                                   insertseq_reads.r1, insertseq_reads.r2,
                                                   c1, c2)

            gr1_mapped = not gr1.is_unmapped
            gr2_mapped = not gr2.is_unmapped
            ir1_mapped = not ir1.is_unmapped
            ir2_mapped = not ir2.is_unmapped


            if gr1_mapped and (not ir1_mapped) and (not gr2_mapped) and ir2_mapped:
                print('G - I')
                gr1_clipped = is_soft_clipped_once(gr1)
                ir2_clipped = is_soft_clipped_once(ir2)
                if (not gr1_clipped) and (not ir2_clipped):

                    gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'G_I', sample)
                    write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                            ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)
                else:
                    pass

            elif (not gr1_mapped) and ir1_mapped and gr2_mapped and (not ir2_mapped):
                print('I - G')
                gr2_clipped = is_soft_clipped_once(gr2)
                ir1_clipped = is_soft_clipped_once(ir1)
                if (not gr2_clipped) and (not ir1_clipped):

                    gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'I_G', sample)
                    write_to_all_flank_sams(gr2, genome_flanks_sam, gr1, genome_flanks_mate_sam,
                                            ir1, insertseq_flanks_sam, ir2, insertseq_flanks_mate_sam)
                    
                else:
                    pass

            elif gr1_mapped and gr2_mapped and (not ir1_mapped) and ir2_mapped:
                gr_concordant = not gr1.mate_is_unmapped
                gr1_clipped = is_soft_clipped_once(gr1)
                gr2_clipped = is_soft_clipped_once(gr2)
                ir2_clipped = is_soft_clipped_once(ir2)

                if gr_concordant:

                    if (not gr1_clipped) and (not gr2_clipped) and (not ir2_clipped):
                        print('G - GI')
                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'G_GI', sample)
                        write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                                ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)

                    elif (not gr1_clipped) and gr2_clipped and ir2_clipped:
                        print("G - GcIc")

                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'G_GcIc', sample)
                        write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                                ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)

                    elif (not gr1_clipped) and (not gr2_clipped) and ir2_clipped:
                        print("G - GIc")

                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'G_GIc', sample)
                        write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                                ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)
                        
                    else:
                        pass
                else:
                    pass

            elif gr1_mapped and gr2_mapped and ir1_mapped and (not ir2_mapped):
                gr_concordant = not gr1.mate_is_unmapped
                gr1_clipped = is_soft_clipped_once(gr1)
                gr2_clipped = is_soft_clipped_once(gr2)
                ir1_clipped = is_soft_clipped_once(ir1)

                if gr_concordant:

                    if (not gr1_clipped) and (not gr2_clipped) and (not ir1_clipped):
                        print("GI - G")

                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'GI_G', sample)
                        write_to_all_flank_sams(gr2, genome_flanks_sam, gr1, genome_flanks_mate_sam,
                                                ir1, insertseq_flanks_sam, ir2, insertseq_flanks_mate_sam)
                        
                    elif gr1_clipped and (not gr2_clipped) and ir1_clipped:
                        print("Gc+Ic - G")

                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'GcIc_G', sample)
                        write_to_all_flank_sams(gr2, genome_flanks_sam, gr1, genome_flanks_mate_sam,
                                                ir1, insertseq_flanks_sam, ir2, insertseq_flanks_mate_sam)
                        
                    elif (not gr1_clipped) and (not gr2_clipped) and ir1_clipped:
                        print("GIc - G")

                        gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'GIc_G', sample)
                        write_to_all_flank_sams(gr2, genome_flanks_sam, gr1, genome_flanks_mate_sam,
                                                ir1, insertseq_flanks_sam, ir2, insertseq_flanks_mate_sam)

                    else:
                        pass
                else:
                    pass

            elif gr1_mapped and ir1_mapped and gr2_mapped and ir2_mapped:

                gr1_clipped = is_soft_clipped_once(gr1)
                gr2_clipped = is_soft_clipped_once(gr2)
                ir1_clipped = is_soft_clipped_once(ir1)
                ir2_clipped = is_soft_clipped_once(ir2)

                ir_concordant = not ir1.mate_is_unmapped

                if ir_concordant:
                    gr_concordant = not gr1.mate_is_unmapped
                    if gr_concordant:
                        if gr1_clipped and gr2_clipped and ir1_clipped and ir2_clipped:

                            if gr1_clipped == gr2_clipped and ir1_clipped == ir2_clipped:
                                print("GcIc - GcIc")

                                gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'GcIc_GcIc', sample)
                                write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                                        ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)

                        elif (not gr1_clipped) and (not gr2_clipped) and ir1_clipped and ir2_clipped:
                            if ir1_clipped == ir2_clipped:
                                print("GIc - GIc")
                                gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'GIc_GIc', sample)
                                write_to_all_flank_sams(gr1, genome_flanks_sam, gr2, genome_flanks_mate_sam,
                                                        ir2, insertseq_flanks_sam, ir1, insertseq_flanks_mate_sam)
            else:
                gr_concordant = not gr1.mate_is_unmapped
                ir_concordant = not ir1.mate_is_unmapped

                gr1, gr2, ir1, ir2 = set_flank_tags(gr1, gr2, ir1, ir2, 'None', sample)

                if gr_concordant:
                    genome_noflanks_sam.write(gr1)
                    genome_noflanks_sam.write(gr2)
                if ir_concordant:
                    insertseq_noflanks_sam.write(ir1)
                    insertseq_noflanks_sam.write(ir2)

            genome_reads = next(genome_sam)
            insertseq_reads = next(insertseq_sam)

        elif c1.name == genome_reads.r1.query_name:
            gr1, gr2 = add_taxon(genome_reads.r1, genome_reads.r2, c1, c2)

            gr_concordant = not gr1.mate_is_unmapped
            if gr_concordant:

                gr1.set_tag('SA', sample)
                gr2.set_tag('SA', sample)

                genome_noflanks_sam.write(gr1)
                genome_noflanks_sam.write(gr2)

            genome_reads = next(genome_sam)


        elif c1.name == insertseq_reads.r1.query_name:
            ir1, ir2 = add_taxon(insertseq_reads.r1, insertseq_reads.r2, c1, c2)

            ir_concordant = not ir1.mate_is_unmapped
            if ir_concordant:

                ir1.set_tag('SA', sample)
                ir2.set_tag('SA', sample)

                insertseq_noflanks_sam.write(ir1)
                insertseq_noflanks_sam.write(ir2)

            insertseq_reads = next(insertseq_sam)

    print("Reads Processed:", count)

if __name__ == '__main__':
    main()