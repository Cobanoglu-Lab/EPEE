import argparse
import math
import os
from time import strftime

# the file generates sbatch script and runs the script over multiple nodes

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--conditiona", help="RNA-seq data for Condition A",
                    type=str, required=True)
parser.add_argument("-b", "--conditionb", help="RNA-seq data for Condition B",
                    type=str, required=True)
parser.add_argument("-na", "--networka", help="Network for condition A",
                    type=str, required=True)
parser.add_argument("-nb", "--networkb", help="Network for condition B",
                    type=str, required=True)
# DEFAULTS
parser.add_argument("-o", "--output", help="output directory", type=str,
                    default='{}_epee_output'.format(strftime('%Y%m%d')))
parser.add_argument("-reg1", "--lregularization", help="lasso regularization \
                    parameter", type=float, default=0.01)
parser.add_argument("-reg2", "--gregularization", help="graph contrained \
                    regularization parameter", type=float, default=0.01)
parser.add_argument("-s", "--step", help="optimizer learning-rate",
                    type=float, default=0.0001)
parser.add_argument("-c", "--conditioning", help="Weight for the interactions \
                    not known", type=bool, default=True)
parser.add_argument("-r", "--runs", help="Number of indpendent runs", type=int,
                    default=14)
parser.add_argument("-i", "--iterations", help="Number of iterations",
                    type=int, default=100000)
parser.add_argument("-norm", "--normalize", help="Normalize the weights",
                    type=str, default='minmax')
parser.add_argument("-model", "--model", help="model to run",
                    type=str, default='epee-gcl')
parser.add_argument("-verbose", "--verbose",
                    help="logging info levels 10, 20, or 30",
                    type=int, default=20)
# OPTIONAL SETTINGS
parser.add_argument("-eval", "--evaluate",
                    help="Evaluation mode available for Th1, Th2, Th17, \
                    Bmem, and AML",
                    type=str, default=None)
parser.add_argument("-prefix", "--prefix",
                    help="Add prefix to the log",
                    type=str, default='')
# OPTIONAL FLAGS
parser.add_argument("-weight", "--weight",
                    help="Store all the inferred weights",
                    action='store_true')
parser.add_argument("-multiprocess", "--multiprocess",
                    help="multiprocess the calculation of perturb and \
                    regulator scores", action='store_true')

# NULL FLAG
parser.add_argument("-null", "--null",
                    help="Generate null scores by shuffling the labels",
                    action='store_true')
# NULL SETTINGS
parser.add_argument("-seed", "--seed", help="Starting seed number",
                    type=int, default=0)
parser.add_argument("-perturb", "--perturb", help="True label perturb scores",
                    type=str, default=None)
# PARAMSCAN SETTINGS
parser.add_argument("-partition", "--partition",
                    help="Partition to run the paramscan",
                    type=str, default='GPU')
parser.add_argument("-nodenumbers", "--nodenumbers",
                    help="Number of nodes",
                    type=int, default=1)
parser.add_argument("-run", "--run",
                    help="Run the sbatch scripts",
                    action='store_true')


if __name__ == '__main__':

    args = parser.parse_args()
    if args.partition == 'GPUv1':
        gres = 'gpu:2'
    else:
        gres = 'gpu:1'

    date = strftime('%Y%m%d')

    # split the params into different nodes
    l1_params = list(map(str, [(1/(10**i)) for i in list(range(0, 9, 1))]))
    l2_range = list(map(str, [(1/(10**i)) for i in list(range(0, 9, 1))]))
    nodes = math.ceil(len(l1_params)/args.nodenumbers)
    l1_params_split = [l1_params[i:i + nodes]
                       for i in range(0, len(l1_params),
                       nodes)]

    print('Number of nodes used: {}'.format(len(l1_params_split)))
    for node, params in enumerate(l1_params_split):
        l1_range = params
        f = open('../sruns/{date}_{context}_paramscan_{n}.sh'
                 .format(date=date, context=args.evaluate, n=node+1), 'w')
        f.write('#!/bin/bash\n')
        f.write('\n')
        f.write('#SBATCH --job-name=scan_{context}_{n}\n'
                .format(context=args.evaluate, n=node+1))

        f.write('#SBATCH --partition={}\n'.format(args.partition))

        f.write('#SBATCH --time=03-00:00:00\n')

        f.write('#SBATCH --output=../../output/log/{date}_{context}_%j.out\n'
                .format(date=date, context=args.evaluate))
        f.write('#SBATCH --error=../../output/log/{date}_{context}_%j.err\n'
                .format(date=date, context=args.evaluate))
        f.write('#SBATCH --mail-type=ALL\n')
        f.write('#SBATCH --mail-user=viren.amin@utsouthwestern.edu\n')
        f.write('\n')
        f.write("cd /project/bioinformatics/Cobanoglu_lab/vamin1/gitlab/epee/"
                "script/sruns\n")
        f.write('\n')

        l1_scan = ' '.join(l1_range)
        l2_scan = ' '.join(l2_range)
        f.write('for i in {}; do\n'.format(l1_scan))
        f.write('\techo $i;\n')
        f.write('\tfor j in {}; do\n'.format(l2_scan))
        f.write('\t\techo $j;\n')
        if args.multiprocess:
            f.write("\t\tCUDA_VISIBLE_DEVICES=0 python ../run_epee_v0.1.4.3.py -a {a} -b {b} -na {na} -nb {nb} -o {o} -reg1 $i -reg2 $j -r {r} -i {i} -model {m} -eval {p} -seed {s} -prefix {prefix} -multiprocess;\n"
                    .format(a=args.conditiona, b=args.conditionb,
                            na=args.networka,
                            nb=args.networkb, o=args.output, r=args.runs,
                            l=args.lregularization, g=args.gregularization,
                            p=args.evaluate, i=args.iterations, m=args.model,
                            s=args.seed, prefix=args.prefix))
        else:
            f.write("\t\tCUDA_VISIBLE_DEVICES=0 python ../run_epee_v0.1.4.3.py -a {a} -b {b} -na {na} -nb {nb} -o {o} -reg1 $i -reg2 $j -r {r} -i {i} -model {m} -eval {p} -seed {s} -prefix {prefix};\n"
                    .format(a=args.conditiona, b=args.conditionb,
                            na=args.networka,
                            nb=args.networkb, o=args.output, r=args.runs,
                            l=args.lregularization, g=args.gregularization,
                            p=args.evaluate, i=args.iterations, m=args.model,
                            s=args.seed, prefix=args.prefix))
        f.write('\tdone;\n')
        f.write('done\n')
        f.close()

        if args.run:
            os.system('sbatch ../sruns/{date}_{context}_paramscan_{n}.sh'
                      .format(date=date, context=args.evaluate, n=node+1))
