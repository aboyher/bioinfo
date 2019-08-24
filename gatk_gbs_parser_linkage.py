import sys, os
import argparse
import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency, fisher_exact


def split_gt_fields(gt):
    return [x.split(":") for x in gt]


def gtfield_sums(gtf, sums):
    for i in range(len(gtf)):
        gt = gtf[i][0]
        if gtf[i] == '.' or gtf[i][0] == './.':
            continue
        if gtf[i][0][0] == '.' and gtf[i][0][0] != gtf[i][0][2]:
            continue
        else:
            sums['total'] += 1
            if gtf[i][0] == '0/0':
                sums['hom0'] += 1
            elif gtf[i][0] == '0/1':
                sums['alt'] += 1
                sums['het'] += 1
            elif gtf[i][0] == '1/1':
                # sums['alt'] += 2
                sums['hom1'] += 1
            # else:
            #     print(gtf[i][0])


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


def parse_variants(fh, inds, qual, cross, outfile, min_parent_prop_het=0.9, sigcutoff=0.001,
                   cmd1=False):
    i = 0
    R1 = np.where([i.startswith("TME14") for i in inds])[0]
    R2 = np.where([i.startswith("NASE14") for i in inds])[0]
    M = np.where([i.startswith("TME204") for i in inds])[0]
    T1 = np.where([i.startswith("5001-26") for i in inds])[0]
    T2 = np.where([i.startswith("5001-46") for i in inds])[0]

    cross_res = cross[cross['CMD (S/R)'] == 'R']
    cross_sus = cross[cross['CMD (S/R)'] == 'S']

    res_idx = np.where([i in cross_res['FullSampleName'].tolist() for i in inds])[0]
    sus_idx = np.where([i in cross_sus['FullSampleName'].tolist() for i in inds])[0]

    outf = open(outfile, 'w')
    print("Res count: {}, Sus count: {}".format(len(res_idx), len(sus_idx)), file = outf)
    # print("Res count: {}, Sus count: {}".format(len(res_idx), len(sus_idx)))
    print("\t\t\t\t\t\t{}".format("\t".join(
        ['R1'] * len(R1) +
        ['R2'] * len(R2) +
        ['T1'] * len(T1) +
        ['T2'] * len(T2) +
        ['r'] * len(res_idx) +
        ['s'] * len(sus_idx))), file = outf)

    # print("\t\t\t\t\t\t{}".format("\t".join(
    #     ['R1'] * len(R1) +
    #     ['R2'] * len(R2) +
    #     ['T1'] * len(T1) +
    #     ['T2'] * len(T2) +
    #     ['r'] * len(res_idx) +
    #     ['s'] * len(sus_idx))))

    print("Contig\tPosition\tmat\tfish_p\thom_mat\thom_fish_p\t{}".format("\t".join(
        np.array(inds)[R1].tolist() +
        np.array(inds)[R2].tolist() +
        np.array(inds)[T1].tolist() +
        np.array(inds)[T2].tolist() +
        np.array(inds)[res_idx].tolist() +
        np.array(inds)[sus_idx].tolist())), file = outf)

    # print("Contig\tPosition\tmat\tfish_p\thom_mat\thom_fish_p\t{}".format("\t".join(
    #     np.array(inds)[R1].tolist() +
    #     np.array(inds)[R2].tolist() +
    #     np.array(inds)[T1].tolist() +
    #     np.array(inds)[T2].tolist() +
    #     np.array(inds)[res_idx].tolist() +
    #     np.array(inds)[sus_idx].tolist())))

    for line in fh:
        i += 1
        variant = line.split()
        # print(variant)

        if float(variant[5]) > qual:
            varout = variant[:2]  # variant[:6]
            varout[0] = "{:12s}".format(varout[0])
            varout[1] = "{:10s}".format(varout[1])
            gtfields = split_gt_fields(variant[9:])

            R1_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[R1]])
            R2_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[R2]])
            T1_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T1]])
            T2_mean = np.mean([x[0] == '0/1' for x in np.array(gtfields)[T2]])
            if not cmd1 and (R1_mean < min_parent_prop_het or T1_mean <
                    min_parent_prop_het or T2_mean < min_parent_prop_het):
                continue
            if cmd1 and (R2_mean <= min_parent_prop_het):
                continue
            #elif not cmd1 and (R1_mean <= min_parent_prop_het):
            #    continue
            #if (T1_mean <= min_parent_prop_het or T2_mean <= min_parent_prop_het):
            #    continue

            gtfields_res = np.array(gtfields)[res_idx].tolist()
            gtfields_sus = np.array(gtfields)[sus_idx].tolist()
            gtfields_ressus = gtfields_res + gtfields_sus

            ## PARSE AND TEST POSITIONS
            if True:
                ## parse each position and test
                res_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}
                sus_sums = {'total': 0, 'alt': 0, 'hom0': 0, 'het': 0, 'hom1': 0}

                gtfield_sums(gtfields_res, res_sums)
                gtfield_sums(gtfields_sus, sus_sums)

                ## TEST EACH ROW HERE
                ## *2 because we are counting both alleles
                # sq_mat = [[res_sums['total'] * 2 - res_sums['alt'], res_sums['alt']],
                #           [sus_sums['total'] * 2 - sus_sums['alt'], sus_sums['alt']]]
                sq_mat = [[res_sums['total'] - res_sums['het'], res_sums['het']],
                          [sus_sums['total'] - sus_sums['het'], sus_sums['het']]]
                # he_mat = [[res_sums['hom0'], res_sums['het']],
                #           [sus_sums['hom']]]
                hom_mat = [[res_sums['hom0'], res_sums['hom1']],
                           [sus_sums['hom0'], sus_sums['hom1']]]

                # print(sq_mat)
                # _, fish_p = fisher_exact(sq_mat)
                fish_p = 0
                _, homfish_p = fisher_exact(hom_mat)
                # _, chi2_p, _, _ = chi2_contingency(sq_mat)
                # theta = (res_sums['alt'] + 2*res_sums['hom1']) / ((res_sums['total'] + sus_sums['total']) * 2)
                # LOD = np.log10( ((1-theta) ** (2*sus_sums['total']) * theta ** (2*res_sums['total']) ) / 0.5 ** ((res_sums['total'] + sus_sums['total']) * 2) )

                # varout.extend(['fishers_p:{:.2e}, hom_fisher_p:{:.2e}, LOD:{:.2e}, mat:{}, hommat:{}'.format(
                #                 fish_p, homfish_p, LOD, sq_mat, hom_mat)])
                varout.extend(['{}\t{:.2e}\t{}\t{:.2e}'.format(sq_mat, fish_p, hom_mat, homfish_p)])

                # if np.argmax([res_sums['hom0'], res_sums['het'], res_sums['hom1']]) != \
                #         np.argmax([sus_sums['hom0'], sus_sums['het'], sus_sums['hom1']]) and \
                #         res_sums['total'] + sus_sums['total'] > 0.5 * len(gtfields_ressus):
                # varout.extend([';'.join(i) for i in gtfields_res])
                # varout.extend([';'.join(i) for i in gtfields_sus])
                # print( "\t".join(varout) )

                # if fish_p <= sigcutoff:
                #     min_count = 4
                #     badvars = []
                #     ## check if there are < min_count resistant 0/0
                #     if hom_mat[0][0] < min_count and hom_mat[0][0] > 0:
                #         gt = np.array([i[0] for i in gtfields_res])
                #         hits = res_idx[gt == '0/0']
                #         badvars.extend( np.array(inds)[hits].tolist() )
                #     if hom_mat[0][1] < min_count and hom_mat[0][1] > 0:
                #         gt = np.array([i[0] for i in gtfields_res])
                #         hits = res_idx[gt == '1/1']
                #         badvars.extend( np.array(inds)[hits].tolist() )
                #     if hom_mat[1][0] < min_count and hom_mat[1][0] > 0:
                #         gt = np.array([i[0] for i in gtfields_sus])
                #         hits = sus_idx[gt == '0/0']
                #         badvars.extend( np.array(inds)[hits].tolist() )
                #     if hom_mat[1][1] < min_count and hom_mat[1][1] > 0:
                #         gt = np.array([i[0] for i in gtfields_sus])
                #         hits = sus_idx[gt == '1/1']
                #         badvars.extend( np.array(inds)[hits].tolist() )
                #
                #     varout.extend([",".join(badvars)])
                if homfish_p <= sigcutoff:
                    ## OUTPUT RES, HET, SUS for each individual and each position
                    for gt in gtfields_res:
                        get_proportion_molecular_correct(varout, res_sums, gt)
                    for gt in gtfields_sus:
                        get_proportion_molecular_correct(varout, sus_sums, gt)

                    print("\t".join(varout), file = outf)
                    # print("\t".join(varout))
                    # sys.exit()
                    # varout.extend(  )

            # print("\t".join(varout))
    print("{} lines".format(i))

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
        print("Error: {} not found".format(args.infile))
        sys.exit(1)

    ## OPEN CROSS FILE
    cross = pd.read_table(args.crossfile, sep='\t')

    if args.cmd1flag:
        cross_interest = ((cross['Cross'] == 'R2xT2') | (cross['Cross'] == 'T2xR2') | (cross['Cross'] == 'T1xR2') | (cross['Cross'] == 'R2xT1') )
    else:
        cross_interest = ((cross['Cross'] == 'R1xT2') | (cross['Cross'] == 'T2xR1') | (cross['Cross'] == 'T1xR1') | (cross['Cross'] == 'R1xT1') | (cross['Cross'] == 'MxT1') | (cross['Cross'] == 'T2xM'))
    cross_t1t2 = cross[cross_interest][['FullSampleName', 'CMD (S/R)']]
    parents = cross[cross['Cross'] == "parental control"]

    outfile = args.infile.split(".vcf")[0] + "_q.{}_p.{}_m.{}_processed.txt".format(args.qual, args.pvalue,
                                                                                    args.min_parent)

    ## OPEN FILE
    with open(args.infile) as fh:
        ## READ COMMENTS IN VCF
        comments, inds = parse_comments(fh)

        ## READ VARIANTS
        parse_variants(fh, inds, float(args.qual), cross_t1t2,outfile,
                       sigcutoff=args.pvalue, min_parent_prop_het=args.min_parent,
                       cmd1=args.cmd1flag)


if __name__ == '__main__':
    main(sys.argv[0], sys.argv[1:])
