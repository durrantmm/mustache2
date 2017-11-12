import sys
import pysam
from numpy import std, mean, median

def read_sam_pairs(sam_path):


    sam_file = pysam.AlignmentFile(sam_path, 'r')

    while True:
        try:
            out = type("", (), {})()
            out.r1 = next(sam_file)
            out.r2 = next(sam_file)
            yield out
        except StopIteration:
            break

def main():
    sam_path = sys.argv[1]
    metrics_path = sys.argv[2]
    counts_path = sys.argv[3]
    sample = sys.argv[4]

    print("Getting insert size stats for %s" % sample)
    sam = read_sam_pairs(sam_path)

    counts = [0]*1001
    tlen_list = list()

    for reads in sam:

        if abs(reads.r1.tlen) == 0 or abs(reads.r1.tlen) > 1000:
            counts[0] += 1
        else:
            tlen_list.append(abs(reads.r1.tlen))
            counts[abs(reads.r1.tlen)] += 1

    with open(metrics_path, 'w') as out:
        out.write('\t'.join([sample, 'MEAN', str(mean(tlen_list))])+'\n')
        out.write('\t'.join([sample, 'MEDIAN', str(median(tlen_list))]) + '\n')
        out.write('\t'.join([sample, 'STDEV', str(std(tlen_list))]) + '\n')

    with open(counts_path, 'w') as out:
        for i in range(len(counts)):
            out.write('\t'.join([sample, str(i), str(counts[i])]))



if __name__ == '__main__':
    main()