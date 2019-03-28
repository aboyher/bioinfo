import sys, os
import argparse
import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency, fisher_exact


def split_gt_fields(gt):
    return [x.split(":") for x in gt]


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


def parse_variants(fh, inds, qual, cross, outfile, mother, min_parent_prop_het=0.9, sigcutoff=0.001):
    inds = [i.split('.R1')[0] for i in inds]
    inds = [i.split('.2')[0] for i in inds]
    inds = [i.split('.3')[0] for i in inds]

    T14 = np.where([i.startswith("TME_14") for i in inds])[0]
    NASE = np.where([i.startswith("NASE_14") for i in inds])[0]
    T04 = np.where([i.startswith("T04026") for i in inds])[0]

    cross_res = cross[cross['Phenotype'] == 'R']
    cross_sus = cross[cross['Phenotype'] == 'S']

    res_idx = np.where([i in cross_res['Sample'].tolist() for i in inds])[0]
    sus_idx = np.where([i in cross_sus['Sample'].tolist() for i in inds])[0]

    outf = open(outfile, 'w')
    print("Res count: {}, Sus count: {}".format(len(res_idx), len(sus_idx)))

    # print("Res count: {}, Sus count: {}".format(len(res_idx), len(sus_idx)))
    print("\t\t\t\t{}".format("\t".join(
        ['r'] * len(res_idx) +
        ['s'] * len(sus_idx))))

    print("Contig\tPosition\thom_mat\thom_fish_p\t{}".format("\t".join(
        np.array(inds)[res_idx].tolist() +
        np.array(inds)[sus_idx].tolist())))

    for line in fh:
        variant = line.split()
        # print(variant)

        if float(variant[5]) > qual:
            varout = variant[:2]  # variant[:6]
            varout[0] = "{:12s}".format(varout[0])
            varout[1] = "{:10s}".format(varout[1])
            gtfields = split_gt_fields(variant[9:])

            NASE_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[NASE]])
            T14_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T14]])
            T04_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T04]])



            # if T04_mean <= min_parent_prop_het or NASE_mean <= min_parent_prop_het or T14_mean <= min_parent_prop_het:
            # continue

            if mother.upper() == "NASE14" and (NASE_mean < 1 or T04_mean < 1):
                continue
            if mother.upper() == "TME14" and (T14_mean < 1 or T04_mean < 1):
                continue
            if mother.upper() == "BOTH" and (T14_mean < 1 or T04_mean < 1 or NASE_mean < 1):
                continue

            gtfields_res = np.array(gtfields)[res_idx].tolist()
            gtfields_sus = np.array(gtfields)[sus_idx].tolist()

            ## PARSE AND TEST POSITIONS
            if True:
                ## parse each position and test
                res_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
                sus_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}

                gtfield_sums(gtfields_res, res_sums)
                gtfield_sums(gtfields_sus, sus_sums)

                hom_mat = [[res_sums['hom0'], res_sums['hom1']],
                           [sus_sums['hom0'], sus_sums['hom1']]]

                fish_p = 0
                _, homfish_p = fisher_exact(hom_mat)

                varout.extend(['{}\t{:.2e}'.format(hom_mat, homfish_p)])

                if homfish_p <= sigcutoff:
                    ## OUTPUT RES, HET, SUS for each individual and each position
                    for gt in gtfields_res:
                        get_proportion_molecular_correct(varout, res_sums, gt)
                    for gt in gtfields_sus:
                        get_proportion_molecular_correct(varout, sus_sums, gt)

                    print("\t".join(varout))
                    # print("\t".join(varout))
                    # sys.exit()
                    # varout.extend(  )

            # print("\t".join(varout))


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


######   MAIN   ######
def main(prog_name, argv):
    # ARG PROCESSING
    parser = argparse.ArgumentParser(prog=prog_name,
                                     description='parse GBS VCF data from CMD2 F1 population',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='infile', type=str, required=True, help='input vcf file')
    parser.add_argument('-c', dest='infofile', type=str, required=True, help='file containing progeny information')
    parser.add_argument('-q', dest='qual', type=int, default=500, help='QUAL score to filter on')
    parser.add_argument('-m', dest='min_parent', type=float, default=0.9,
                        help='min proportion of parents required to be het')
    parser.add_argument('-p', dest='pvalue', type=float, default=0.001,
                        help='limit output to max pvalue from fisher\'s exact test')
    parser.add_argument('-M', dest='mother', type=str, default='BOTH')
    args = parser.parse_args(argv)

    ## check for file
    if not os.path.isfile(args.infile):
        print("Error: {} not found".format(args.infile))
        sys.exit(1)

    ## OPEN CROSS FILE

    prog = pd.read_table(args.infofile, sep='\t')

    if args.mother.upper() == 'BOTH':
        prog_interest = ((prog['Mother'] == "NASE 14") | (prog['Mother'] == "TME 14"))
    elif args.mother.upper() == 'NASE14':
        prog_interest = ((prog["Mother"] == "NASE 14"))
    elif args.mother.upper() == 'TME14':
        prog_interest = ((prog["Mother"] == "TME 14"))
    else:
        print("Error: {} not a recognized mother. Options are 'both', 'NASE14', 'TME14'")
    prog_interest = prog[prog_interest][["Sample", "Phenotype"]]


    #
    outfile = args.infile.split(".vcf")[0] + "_q.{}_p.{}_m.{}_" + args.mother + "_processed.txt".format(args.qual, args.pvalue,
                                                                                    args.min_parent)
    #
    # ## OPEN FILE
    with open(args.infile) as fh:
        ## READ COMMENTS IN VCF
        comments, inds = parse_comments(fh)

        ## READ VARIANTS
        parse_variants(fh, inds, float(args.qual), prog_interest, outfile, args.mother,
                       sigcutoff=args.pvalue, min_parent_prop_het=args.min_parent)


if __name__ == '__main__':
    main(sys.argv[0], sys.argv[1:])
