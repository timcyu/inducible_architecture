import sys
import csv
import math
import numpy as np

def calculateUninducedExpression(E10, E35, Eprox, Edist, Eloop, r_min = 0.879, r_max = 25.14, e = 2.71828):
	# calculate the boltzmann weights
    Ecore = E10 + E35
    prob_core_prom = e**-Ecore
    prob_dist = e**-Edist
    prob_prox = e**-Eprox
    prob_prox_dist = e**-(Eprox + Edist)
    prob_prox_dist_loop = e**-(Eprox + Edist + Eloop)

    # formula for gene expression
    theta = (1 + prob_core_prom)*(1 + prob_dist + prob_prox + prob_prox_dist + prob_prox_dist_loop)
    gene_exp = (r_min*(theta - (prob_core_prom*(1 + prob_dist))) + r_max*(prob_core_prom*(1 + prob_dist)))/theta

    return gene_exp


def calculateInducedExpression(E10, E35, Eprox, Edist, Eloop, r_min = 0.879, r_max = 25.14, e = 2.71828, p_actrep = 0.0282):
	# these relationships allow us to model expression at 1mM IPTG
    ind_Eprox = Eprox - math.log(p_actrep)
    ind_Edist = Edist - math.log(p_actrep)
    ind_Eloop = Eloop + math.log(p_actrep)

    # calculate the boltzmann weights
    Ecore = E10 + E35
    prob_core_prom = e**-Ecore
    prob_dist = e**-ind_Edist
    prob_prox = e**-ind_Eprox
    prob_prox_dist = e**-(ind_Eprox + ind_Edist)
    prob_prox_dist_loop = e**-(ind_Eprox + ind_Edist + ind_Eloop)

    # formula for gene expressiion
    theta = (1 + prob_core_prom)*(1 + prob_dist + prob_prox + prob_prox_dist + prob_prox_dist_loop)
    gene_exp = (r_min*(theta - (prob_core_prom*(1 + prob_dist))) + r_max*(prob_core_prom*(1 + prob_dist)))/theta

    return gene_exp

if __name__ == "__main__":

    # initialize parameters
    Eloop = -2.49                
    #E10s = [0, 4, 5.03, 3.77] 
    #E35s = [-1.96, 3.58, 0.66, 1.7]
    #Eproxs = [-0.23, 3.01, 3.91, 27.88, 1.58, 2.01, 3.76, 25.96, -2.09]
    #Edists = [-0.23, 3.01, 3.91, 27.88, 1.58, 2.01, 3.76, 25.96, -2.09]

    E10s = [0, 3.77] # consensus -10 and strongest mutated 
    E35s = [-1.96, 0.66] # consensus -35 and strongest mutated
    Eproxs = np.linspace(-2.2, 8, num = 100, endpoint = False).tolist()
    Edists = [-0.23, 3.01, 3.91, 27.88, 1.58, 2.01, 3.76, 25.96, -2.09]

    # open an output file to store the data
    file = open("modelExp.txt", "w")
    file.write("induced_exp" + '\t' + "uninduced_exp" + '\t' + "fold_change" + '\t' + "E10" + '\t' + "E35" + '\t' + "Eprox" + '\t' + "Edist" + '\t' + "Condition" + '\n')

    for i in E10s:
        for j in E35s:
            for k in Eproxs:
                for l in Edists:
                    unind_exp = calculateUninducedExpression(i, j, k, l, Eloop)
                    ind_exp = calculateInducedExpression(i, j, k, l, Eloop)
                    fold_change = ind_exp/unind_exp
                    file.write(str(ind_exp) + '\t' + str(unind_exp) + '\t' + str(fold_change) + '\t' + str(i) + '\t' + str(j) + '\t' + str(k) + '\t' + str(l) + '\t' + 'FC' + '\n')
    file.close()




