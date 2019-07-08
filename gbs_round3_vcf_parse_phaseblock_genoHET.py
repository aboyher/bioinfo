import sys, os
import argparse
import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency, fisher_exact

def parse_comments(fh):
    comments = []
    inds = []
    line = fh.readline()

    # sys.stdout.write("start comments\n")
    while line.startswith("##"):
        line = fh.readline()
        # sys.stdout.write("comment\n")
        # pass
        comments.append(line)
    # sys.stdout.write("done comments\n")

    inds = line.split()[9:]

    return comments, inds

def read_cross(fh):
    cross = pd.read_table()
    cross = {}

    ## skip first line
    fh.readline()

    for line in fh:
        line = line.split()
        cross[line[6]] = [line[10]] + line[7].split('x')

    return cross

def get_proportion_molecular_correct(varout, sums, gt):
    if gt[0] == '0/1':
        varout.append("het")
    else:
        try:
            homs = [sums['hom0'], sums['hom1']]
            if gt[0] == '0/0':
                varout.append(str(homs[0] / float(sum(homs))))
            elif gt[0] == '1/1':
                varout.append(str(homs[1] / float(sum(homs))))
            else:
                varout.append("N/A")
        except:
            varout.append("N/A")

def get_genotype(varout, sums, gt):
    try:
        gts = gt[0].split("/")
        if gts[0] == gts[1]:
            varout.append(gts[0])
        else:
            varout.append("het")
    except:
        varout.append("n/a")

def gtfield_sums(gtf, sums):
    for i in range(len(gtf)):
        if gtf[i] == '.' or gtf[i][0] == '.':
            continue
        else:
            sums['total'] += 1
            if gtf[i][0] == '0/0':
                sums['hom0'] += 1
            elif gtf[i][0] == '0/1':
                sums['alt'] += 1
                sums['het'] += 1
            elif gtf[i][0] == '1/1':
                sums['hom1'] += 1


def split_gt_fields(gt):
    return [x.split(":") for x in gt]

def parse_variants(fh, inds , cross, outfile, cmd1=False):
    inds = [i.split('.R1')[0] for i in inds]
    T1 = np.where([i.startswith("5001-26") for i in inds])[0]
    T2 = np.where([i.startswith("5001-46") for i in inds])[0]
    R1 = np.where([i.startswith("TME14") for i in inds])[0]
    R2 = np.where([i.startswith("NASE14") for i in inds])[0]

    cross_res = cross[cross['CMD (S/R)'] == 'R']
    cross_sus = cross[cross['CMD (S/R)'] == 'S']
    res_idx = np.where([i in cross_res['FullSampleName'].tolist() for i in inds])[0]
    sus_idx = np.where([i in cross_sus['FullSampleName'].tolist() for i in inds])[0]

    outf = open(outfile,'w')
    print("\t\t{}".format("\t".join(
    ['R1']*len(R1) +
    ['R2']*len(R2) +
    ['T1']*len(T1) +
    ['T2']*len(T2) +
    ['r']*len(res_idx) +
    ['s']*len(sus_idx))), file = outf)
    print("Scaffold\tPosition\t{}".format("\t".join(
    np.array(inds)[R1].tolist() +
    np.array(inds)[R2].tolist() +
    np.array(inds)[T1].tolist() +
    np.array(inds)[T2].tolist() +
    np.array(inds)[res_idx].tolist() +
    np.array(inds)[sus_idx].tolist())), file = outf)
    #print("R1: {}, R2:{}, T1:{}, T2:{}, Res:{}, Sus:{}\t{}".format(
    #    len(R1),
    #    len(R2),
    #    len(T1),
    #    len(T2),
    #    len(res_idx),
    #    len(sus_idx),
    #    "\t".join(
            # np.array(inds)[R1].tolist() +
            # np.array(inds)[R2].tolist() +
            # np.array(inds)[T1].tolist() +
            # np.array(inds)[T2].tolist() +
            # np.array(inds)[res_idx].tolist() +
            # np.array(inds)[sus_idx].tolist())))

    for line in fh:
        variant = line.split()
       # if variant[0] != 'Super-Scaffold_1' and variant[0] != 'Super-Scaffold_26':
       #     continue
       # if np.float(variant[5]) < 1000:
       #     continue

        varout = variant[:2]
        varout[0] = "{:12s}".format(varout[0])
        varout[1] = "{:10s}".format(varout[1])
        #varout.append("{:10s}".format(variant[7].split(";")[7].split("=")[1]))
        gtfields = split_gt_fields(variant[9:])

        R1_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[R1]])
        R2_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[R2]])
        T1_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T1]])
        T2_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T2]])
        if R1_mean != 1 or T1_mean != 1 or T2_mean != 1:
            continue
        # print ("R1: {} T1: {} T2: {}".format([x[0] for x in np.array(gtfields)[R1]],[x[0] for x in np.array(gtfields)[T1]], [x[0] for x in np.array(gtfields)[T2]]))

        gtfields_R1 = np.array(gtfields)[R1].tolist()
        gtfields_R2 = np.array(gtfields)[R2].tolist()
        gtfields_T1 = np.array(gtfields)[T1].tolist()
        gtfields_T2 = np.array(gtfields)[T2].tolist()
        gtfields_res = np.array(gtfields)[res_idx].tolist()
        gtfields_sus = np.array(gtfields)[sus_idx].tolist()


        R1_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
        R2_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
        T1_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
        T2_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
        res_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
        sus_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}

        gtfield_sums(gtfields_R1, R1_sums)
        gtfield_sums(gtfields_R2, R2_sums)
        gtfield_sums(gtfields_T1, T1_sums)
        gtfield_sums(gtfields_T2, T2_sums)
        gtfield_sums(gtfields_res, res_sums)
        gtfield_sums(gtfields_sus, sus_sums)

        hom_mat = [[res_sums['hom0'], res_sums['hom1']],[sus_sums['hom0'], sus_sums['hom1']]]


        _, fish_p = fisher_exact(hom_mat)
        #varout.extend(['{:2e}\t{}'.format(fish_p, hom_mat)])

        # for gt in gtfields_R1:
        #     get_proportion_molecular_correct(varout, R1_sums, gt)
        # for gt in gtfields_R2:
        #     get_proportion_molecular_correct(varout, R2_sums, gt)
        # for gt in gtfields_T1:
        #     get_proportion_molecular_correct(varout, T1_sums, gt)
        # for gt in gtfields_T2:
        #     get_proportion_molecular_correct(varout, T2_sums, gt)
        # for gt in gtfields_res:
        #     get_proportion_molecular_correct(varout, res_sums, gt)
        # for gt in gtfields_sus:
        #     get_proportion_molecular_correct(varout, sus_sums, gt)


        for gt in gtfields_R1:
            get_genotype(varout, R1_sums, gt)
        for gt in gtfields_R2:
            get_genotype(varout, R2_sums, gt)
        for gt in gtfields_T1:
            get_genotype(varout, T1_sums, gt)
        for gt in gtfields_T2:
            get_genotype(varout, T2_sums, gt)
        for gt in gtfields_res:
            get_genotype(varout, res_sums, gt)
        for gt in gtfields_sus:
            get_genotype(varout, sus_sums, gt)

        print("\t".join(varout), file = outf)


def main(prog_name, argv):
    # ARG PROCESSING
    parser = argparse.ArgumentParser(prog=prog_name,
                                     description='parse GBS VCF data from CMD2 F1 population',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='infile', type=str, required=True, help='input vcf file')
    parser.add_argument('-c', dest='crossfile', type=str, required=True, help='file containing cross information')
    parser.add_argument('-q', dest='qual', type=int, default=500, help='QUAL score to filter on')
    parser.add_argument('-m', dest='min_parent', type=float, default=0.9,
                        help='min proportion of parents required to be het')
    parser.add_argument('-p', dest='pvalue', type=float, default=0.001,
                        help='limit output to max pvalue from fisher\'s exact test')
    parser.add_argument('--cmd1', dest='cmd1flag', action='store_true',
                        help='Flag to use CMD1 types instead of CMD2 types')
    # parser.add_argument('-', nargs='+', type=str, help='mpileup output files to read from')
    args = parser.parse_args(argv)

    ## check for file
    if not os.path.isfile(args.infile):
        sys.exit(1)
        print("Error: {} not found".format(args.infile))

    ## OPEN CROSS FILE
    cross = pd.read_table(args.crossfile, sep='\t')

    if args.cmd1flag:
        cross_interest = ((cross['Cross'] == 'R2xT2') | (cross['Cross'] == 'T2xR2') | (cross['Cross'] == 'T1xR2') | (cross['Cross'] == 'R2xT1') )
    else:
        cross_interest = ((cross['Cross'] == 'R1xT2') | (cross['Cross'] == 'T2xR1') | (cross['Cross'] == 'T1xR1') | (cross['Cross'] == 'R1xT1') | (cross['Cross'] == 'MxT1') | (cross['Cross'] == 'T2xM'))
    cross_t1t2 = cross[cross_interest][['FullSampleName', 'CMD (S/R)']]
    parents = cross[cross['Cross'] == "parental control"]

    outfile = args.infile.split(".vcf")[0] + "_geno_parentHET.txt"

    ## OPEN FILE
    with open(args.infile) as fh:
        ## READ COMMENTS IN VCF
        comments, inds = parse_comments(fh)

        ## READ VARIANTS
        parse_variants(fh, inds, cross_t1t2,outfile, cmd1=args.cmd1flag)


if __name__ == '__main__':
    main(sys.argv[0], sys.argv[1:])
