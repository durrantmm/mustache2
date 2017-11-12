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

def adjust_read_order_once(r1, r2, class1, class2):

    if r1.is_read1:
        r1.set_tag('TX', class1.taxon)
        r2.set_tag('TX', class2.taxon)
        out_gr = [r1, r2]
    else:
        r1.set_tag('TX', class2.taxon)
        r2.set_tag('TX', class1.taxon)
        out_gr = [r2, r1]

    return out_gr

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

def all_reads_match(c1, genome_reads, insertseq_reads):
    return c1.name == genome_reads.r1.query_name == insertseq_reads.r1.query_name


def genome_reads_match(c1, genome_reads):
    return c1.name == genome_reads.r1.query_name


def insertseq_reads_match(c1, insertseq_reads):
    return c1.name == insertseq_reads.r1.query_name

def set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, scenario):
    gr1.set_tag('SC', scenario), gr2.set_tag('SC', scenario), ir1.set_tag('SC', scenario), ir1.set_tag('SC', scenario)
    gr1.set_tag('SP', sample), gr2.set_tag('SP', sample), ir1.set_tag('SP', sample), ir1.set_tag('SP', sample)
    gr1.set_tag('IS', ir2.reference_name), gr2.set_tag('IS', ir1.reference_name)
    ir1.set_tag('GN', gr2.reference_name), ir1.set_tag('GN', gr1.reference_name)

def check_scenario1(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 1")
    if gr1_mapped and (not gr2_mapped) and (not ir1_mapped) and ir2_mapped:
        print("SCENARIO 1A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 1)

        genome_flanks_sam.write(gr1)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    elif (not gr1_mapped) and gr2_mapped and ir1_mapped and (not ir2_mapped):
        print("SCENARIO 1B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 1)

        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)

        return True

    else:
        return False


def check_scenario2(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 2")

    gr_share_contig = gr1.reference_name == gr2.reference_name
    if gr1_mapped and gr2_mapped and (not gr_share_contig):

        if (not ir1_mapped) and ir2_mapped:
            print("SCENARIO 2A")
            set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 2)

            genome_flanks_sam.write(gr1)
            genome_flanks_mate_sam.write(gr2)
            insertseq_flanks_sam.write(ir2)
            insertseq_flanks_mate_sam.write(ir1)

            return True

        elif ir1_mapped and (not ir2_mapped):
            print("SCENARIO 2B")
            set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 2)

            genome_flanks_sam.write(gr2)
            genome_flanks_mate_sam.write(gr1)
            insertseq_flanks_sam.write(ir1)
            insertseq_flanks_mate_sam.write(ir2)

            return True
    else:
        return False


def check_scenario3(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 3")
    gr_distantly_mapped = False


    if abs(gr1.pos - gr2.pos) > 5000:
        #print(gr1.pos, gr2.pos, abs(gr1.pos - gr2.pos))
        gr_distantly_mapped = True

    if gr1_mapped and gr2_mapped and gr_distantly_mapped:

        if (not ir1_mapped) and ir2_mapped:
            print("SCENARIO 3A")
            set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 3)
            genome_flanks_sam.write(gr1)
            genome_flanks_mate_sam.write(gr2)
            insertseq_flanks_sam.write(ir2)
            insertseq_flanks_mate_sam.write(ir1)

            return True

        elif ir1_mapped and (not ir2_mapped):
            print("SCENARIO 3B")
            set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 3)
            genome_flanks_sam.write(gr2)
            genome_flanks_mate_sam.write(gr1)
            insertseq_flanks_sam.write(ir1)
            insertseq_flanks_mate_sam.write(ir2)

            return True
    else:
        return False

def check_scenario4(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 3")

    if (not ir1_mapped) or (not ir2_mapped):
        return False

    ir_share_contig = ir1.reference_name == ir2.reference_name
    if not ir_share_contig:
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring and gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring and gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring.count('S'): ir2_soft_clipped_once = True


    if (gr1_mapped and gr1_soft_clipped_once) and (not gr2_mapped) and (ir1_mapped and ir1_soft_clipped_once) \
            and (ir2_mapped and not ir2_soft_clipped_once):

        if (not ir1.is_reverse) and ir2.is_reverse:
            if ir1.cigartuples[0][0] != 4:
                return False
        elif ir1.is_reverse and (not ir2.is_reverse):
            if ir1.cigartuples[-1][0] != 4:
                return False

        print("SCENARIO 4A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 4)
        genome_flanks_sam.write(gr1)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    elif (not gr1_mapped) and (gr2_mapped and gr2_soft_clipped_once) and (ir1_mapped and not ir1_soft_clipped_once) \
            and (ir2_mapped and ir2_soft_clipped_once):

        if ir1.is_reverse and (not ir2.is_reverse):
            if ir2.cigartuples[0][0] != 4:
                return False
        elif (not ir1.is_reverse) and ir2.is_reverse:
            if ir1.cigartuples[-1][0] != 4:
                return False

        print("SCENARIO 4B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 4)
        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)


        return True

    else:
        return False


def check_scenario5(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 3")

    if (not gr1_mapped) or (not gr2_mapped):
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring and ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring and ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    if (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and gr2_soft_clipped_once) and \
            (not ir1_mapped) and (ir2_mapped and ir2_soft_clipped_once):

        if gr1.is_reverse and (not gr2.is_reverse):
            if gr2.cigartuples[0][0] != 4:
                return False
        elif (not gr1.is_reverse) and gr2.is_reverse:
            if gr1.cigartuples[-1][0] != 4:
                return False

        print("SCENARIO 5A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 5)
        genome_flanks_sam.write(gr1)
        #genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    elif (gr1_mapped and gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) and \
            (ir1_mapped and ir1_soft_clipped_once) and (not ir2_mapped):

        if (not gr1.is_reverse) and gr2.is_reverse:
            if gr1.cigartuples[0][0] != 4:
                return False
        elif gr1.is_reverse and (not gr2.is_reverse):
            if gr1.cigartuples[-1][0] != 4:
                return False


        print("SCENARIO 5B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 5)
        #genome_flanks_sam.write(gr1)
        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)

        return True

    else:
        return False

def check_scenario6(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 6")

    if (not gr1_mapped) or (not gr2_mapped) or (not ir1_mapped) or (not ir2_mapped):
        return False

    ir_share_contig = ir1.reference_name == ir2.reference_name
    if not ir_share_contig:
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    if gr1_soft_clipped_once and gr2_soft_clipped_once and ir1_soft_clipped_once and ir2_soft_clipped_once:

        gr_left_softclipped = gr1.cigartuples[0] == gr2.cigartuples[0]
        gr_right_softclipped = gr1.cigartuples[-1] == gr2.cigartuples[-1]

        ir_left_softclipped = ir1.cigartuples[0] == ir2.cigartuples[0]
        ir_right_softclipped = ir1.cigartuples[-1] == ir2.cigartuples[-1]

        if (not gr_left_softclipped) and (not gr_right_softclipped):
            return False
        if (not ir_left_softclipped) and (not ir_right_softclipped):
            return False

        #print(gr1.cigartuples, gr1.query_sequence)
        #print(gr2.cigartuples, gr2.query_sequence)
        #print(ir1.cigartuples, ir1.query_sequence)
        #print(ir2.cigartuples, ir2.query_sequence)

        if gr1.query_sequence == ir1.query_sequence:
            #print("SAME SEQS")
            if (gr_right_softclipped and ir_right_softclipped) or \
                (gr_left_softclipped and ir_left_softclipped):
                #"SOFTCLIPPED ON THE SAME SIDE"
                return False
        else:
            #print("DIFF SEQS")
            if (gr_right_softclipped and ir_left_softclipped) or \
                    (gr_left_softclipped and ir_right_softclipped):
                # "SOFTCLIPPED ON THE SAME SIDE"
                return False

        print("SCENARIO 6")

        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 6)

        gr_flank, gr_mate, ir_flank, ir_mate = gr1, gr2, gr2, gr1
        if gr_left_softclipped:
            if gr_flank.cigartuples[0][1] > gr_mate.cigartuples[0][1]:
                gr_flank, gr_mate, ir_flank, ir_mate = gr2, gr1, gr1, gr2
        if gr_right_softclipped:
            if gr_flank.cigartuples[-1][1] > gr_mate.cigartuples[-1][1]:
                gr_flank, gr_mate, ir_flank, ir_mate = gr2, gr1, gr1, gr2

        genome_flanks_sam.write(gr_flank)
        #genome_flanks_sam.write(gr_mate)
        genome_flanks_mate_sam.write(gr_mate)
        insertseq_flanks_sam.write(ir_flank)
        insertseq_flanks_mate_sam.write(ir_mate)

        return True

    else:
        return False

def check_scenario7(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 1")

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring and gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring and gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring and ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring and ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    #if (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) and \
    #        (not ir1_mapped) and (ir2_mapped and not ir2_soft_clipped_once):
    if gr1_mapped and gr2_mapped  and not ir1_mapped and ir2_mapped and not ir2_soft_clipped_once:
        print("SCENARIO 7A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 7)

        genome_flanks_sam.write(gr1)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    #elif (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) and \
    #        (ir1_mapped and not ir1_soft_clipped_once) and (not ir2_mapped):
    elif gr1_mapped and gr2_mapped and ir1_mapped and (not ir2_mapped):
        print("SCENARIO 7B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 7)

        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)

        return True

    else:
        return False


def check_scenario8(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 3")

    if (not ir1_mapped) or (not ir2_mapped):
        return False

    ir_share_contig = ir1.reference_name == ir2.reference_name
    if not ir_share_contig:
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring and gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring and gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    if (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) \
            and (ir1_mapped and ir1_soft_clipped_once) and (ir2_mapped and not ir2_soft_clipped_once):

        if (not ir1.is_reverse) and ir2.is_reverse:
            if ir1.cigartuples[0][0] != 4:
                return False
        elif ir1.is_reverse and (not ir2.is_reverse):
            if ir1.cigartuples[-1][0] != 4:
                return False

        print("SCENARIO 8A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 8)
        genome_flanks_sam.write(gr1)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    elif (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and gr2_soft_clipped_once) \
            and (ir1_mapped and not ir1_soft_clipped_once) and (ir2_mapped and ir2_soft_clipped_once):

        if ir1.is_reverse and (not ir2.is_reverse):
            if ir2.cigartuples[0][0] != 8:
                return False
        elif (not ir1.is_reverse) and ir2.is_reverse:
            if ir1.cigartuples[-1][0] != 8:
                return False

        print("SCENARIO 8B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 8)
        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)


        return True

    else:
        return False

def check_scenario9(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 3")

    if (not gr1_mapped) or (not gr2_mapped):
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring and ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring and ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    if (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) and \
            (not ir1_mapped) and (ir2_mapped and ir2_soft_clipped_once):

        if gr1.is_reverse and (not gr2.is_reverse):
            if gr2.cigartuples[0][0] != 4:
                return False
        elif (not gr1.is_reverse) and gr2.is_reverse:
            if gr1.cigartuples[-1][0] != 4:
                return False

        print("SCENARIO 9A")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 9)
        genome_flanks_sam.write(gr1)
        #genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    elif (gr1_mapped and not gr1_soft_clipped_once) and (gr2_mapped and not gr2_soft_clipped_once) and \
            (ir1_mapped and ir1_soft_clipped_once) and (not ir2_mapped):

        if (not gr1.is_reverse) and gr2.is_reverse:
            if gr1.cigartuples[0][0] != 4:
                return False
        elif gr1.is_reverse and (not gr2.is_reverse):
            if gr1.cigartuples[-1][0] != 4:
                return False


        print("SCENARIO 9B")
        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 9)
        #genome_flanks_sam.write(gr1)
        genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr1)
        insertseq_flanks_sam.write(ir1)
        insertseq_flanks_mate_sam.write(ir2)

        return True

    else:
        return False


def check_scenario10(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped,
                    genome_flanks_sam, genome_flanks_mate_sam, insertseq_flanks_sam, insertseq_flanks_mate_sam):
    # SCENARIO 1A
    # READ 1 ALIGNS TO GENOME, READ 2 ALIGNS TO IS
    #print("SCENARIO 6")

    if (not gr1_mapped) or (not gr2_mapped) or (not ir1_mapped) or (not ir2_mapped):
        return False

    ir_share_contig = ir1.reference_name == ir2.reference_name
    if not ir_share_contig:
        return False

    gr1_soft_clipped_once, gr2_soft_clipped_once, ir1_soft_clipped_once, ir2_soft_clipped_once, = False, False, False, False
    if gr1.cigarstring.count('S'): gr1_soft_clipped_once = True
    if gr2.cigarstring.count('S'): gr2_soft_clipped_once = True
    if ir1.cigarstring.count('S'): ir1_soft_clipped_once = True
    if ir2.cigarstring.count('S'): ir2_soft_clipped_once = True

    if (not gr1_soft_clipped_once) and (not gr2_soft_clipped_once) and ir1_soft_clipped_once and ir2_soft_clipped_once:

        gr_left_softclipped = gr1.cigartuples[0] == gr2.cigartuples[0]
        gr_right_softclipped = gr1.cigartuples[-1] == gr2.cigartuples[-1]

        ir_left_softclipped = ir1.cigartuples[0] == ir2.cigartuples[0]
        ir_right_softclipped = ir1.cigartuples[-1] == ir2.cigartuples[-1]

        if (not gr_left_softclipped) and (not gr_right_softclipped):
            return False
        if (not ir_left_softclipped) and (not ir_right_softclipped):
            return False

        print("SCENARIO 10")

        set_final_tags(sample, c1, c2, gr1, gr2, ir1, ir2, 10)

        genome_flanks_sam.write(gr1)
        #genome_flanks_sam.write(gr2)
        genome_flanks_mate_sam.write(gr2)
        insertseq_flanks_sam.write(ir2)
        insertseq_flanks_mate_sam.write(ir1)

        return True

    else:
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

        #print(c1.name, c2.name, genome_reads.r1.query_name)
        # BEGIN PROCESSING READ INFORMATION
        if all_reads_match(c1, genome_reads, insertseq_reads):
            #print("ALL READS MATCH")
            gr1, gr2, ir1, ir2 = adjust_read_order(genome_reads.r1, genome_reads.r2,
                                                   insertseq_reads.r1, insertseq_reads.r2,
                                                   c1, c2)

            gr1_mapped, gr2_mapped = (not gr1.is_unmapped), (not gr2.is_unmapped)
            ir1_mapped, ir2_mapped = (not ir1.is_unmapped), (not ir2.is_unmapped)

            #print(gr1_mapped, gr2_mapped, ir1_mapped, ir2_mapped)
            scenario1_applies = check_scenario1(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario1_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario2_applies = check_scenario2(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario2_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario3_applies = check_scenario3(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario3_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario4_applies = check_scenario4(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)



            if scenario4_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario5_applies = check_scenario5(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario5_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario6_applies = check_scenario6(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario6_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario7_applies = check_scenario7(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)
            if scenario7_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario8_applies = check_scenario8(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario8_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario9_applies = check_scenario9(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                insertseq_flanks_sam, insertseq_flanks_mate_sam)

            if scenario9_applies:
                genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)
                continue

            scenario10_applies = check_scenario10(sample, c1, c2, gr1, gr2, ir1, ir2, gr1_mapped, gr2_mapped, ir1_mapped,
                                                  ir2_mapped, genome_flanks_sam, genome_flanks_mate_sam,
                                                  insertseq_flanks_sam, insertseq_flanks_mate_sam)

            genome_reads, insertseq_reads = next(genome_sam), next(insertseq_sam)


        elif genome_reads_match(c1, genome_reads):
            #print("GENOME READS MATCH")
            gr1, gr2 = adjust_read_order_once(genome_reads.r1, genome_reads.r2, c1, c2)
            gr1.set_tag('SP', sample), gr2.set_tag('SP', sample)

            genome_noflanks_sam.write(gr1)
            genome_noflanks_sam.write(gr2)

            genome_reads = next(genome_sam)

        elif insertseq_reads_match(c1, insertseq_reads):
            #print("INSERTSEQ READS MATCH")
            ir1, ir2 = adjust_read_order_once(insertseq_reads.r1, insertseq_reads.r2, c1, c2)
            ir1.set_tag('SP', sample), ir2.set_tag('SP', sample)

            insertseq_noflanks_sam.write(ir1)
            insertseq_noflanks_sam.write(ir2)

            insertseq_reads = next(insertseq_sam)

    print("Reads Processed:", count)

if __name__ == '__main__':
    main()