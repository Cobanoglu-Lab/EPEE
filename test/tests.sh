# # testing scenario 1 without multiprocessing
# python ../script/run_epee_v0.1.4.3.py --conditiona data/CD4_Naive.txt.gz --conditionb data/CD4_Th2.txt.gz --networka data/cd4+_t_cells.txt.gz --networkb data/cd4+_t_cells.txt.gz -r 2 -i 100 -prefix test_out
#
# # testing scenario 2 with multiprocessing
# python ../script/run_epee_v0.1.4.3.py --conditiona data/CD4_Naive.txt.gz --conditionb data/CD4_Th2.txt.gz --networka data/cd4+_t_cells.txt.gz --networkb data/cd4+_t_cells.txt.gz -r 2 -i 100 -prefix test_out -multiprocess
#
# # testing scenario 3 null without multiprocessing
# python ../script/run_epee_v0.1.4.3.py --conditiona data/CD4_Naive.txt.gz --conditionb data/CD4_Th2.txt.gz --networka data/cd4+_t_cells.txt.gz --networkb data/cd4+_t_cells.txt.gz -r 2 -i 100 -prefix test_out -perturb data/perturb_score.txt.gz -null
#
# # testing scenario 4 null with multiprocessing
# python ../script/run_epee_v0.1.4.3.py --conditiona data/CD4_Naive.txt.gz --conditionb data/CD4_Th2.txt.gz --networka data/cd4+_t_cells.txt.gz --networkb data/cd4+_t_cells.txt.gz -r 2 -i 100 -prefix test_out -perturb data/perturb_score.txt.gz -null -multiprocess

# testing tfwrapper
../script/analysis/tfwrapper_paramsrunner.sh python ../script/run_epee_v0.1.4.3.py --conditiona data/CD4_Naive.txt.gz --conditionb data/CD4_Th2.txt.gz --networka data/cd4+_t_cells.txt.gz --networkb data/cd4+_t_cells.txt.gz -r 2 -i 100 -prefix test_tfwrapper
