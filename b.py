import matplotlib.pyplot as plt
from numpy import*
import numpy as np
#from fsi_file import simInt_fsi
#from geneScore_geneInt import Getscorepergene, GetGeneInteraction
#from generateFile_si import getInfos, generate_file

import pandas as pd
import openpyxl
import os
from openpyxl.utils.dataframe import dataframe_to_rows

import csv
import copy
from threading import Thread
from random import Random
from random import sample
from time import time
import inspyred
from inspyred import ec
from inspyred.ec import terminators
import sqlite3 as lite
import sys
#from Tkinter import *
from operator import itemgetter
import itertools
import random
from datetime import datetime
from multiprocessing import Manager 

import logging
import logging.handlers #Handlers send the log records (created by loggers) to the appropriate destination.
logger = logging.getLogger('inspyred.ec')
logger.propagate = False 
logger.setLevel(logging.DEBUG)
file_handler = logging.handlers.RotatingFileHandler('inspyredetest.log', mode='w')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

def main(display=False):
    
    start_time = datetime.now()
    rand = Random() 
    tabmax = 40 # nbre max de gènes qu'une communauté peut avoir = 40

    def getsizeCommunitie(): 
        var = rand.randint(5,tabmax) # taille de communauté choisie aléatoirement (5,40)
        return var
    
    # créer une liste de gènes à partir de "pathW2.txt"                
    def getGeneFromFile():

        f2 = open ('pathW2.txt','r')
        l = list() 
        
        for  i in  f2.readlines(): 
            lignes = i.strip('\n ') 
            l.append(lignes) # l= liste contenant tous les gènes obtenus du fichier pathW2.txt
        
        nv = list() # communauté
        v = getsizeCommunitie()
        
        # créer la communauté sans redondance de gènes à partir de la liste l.
        for i in range(0,v):
            i = rand.choice(l) 
            l.remove(i) 
            nv.append(i)
        return nv       
    
    # générer pop initiale à partir d'un fichier
    def generate_pop(random, args): 
        pop = getGeneFromFile()
        print (pop)
        return pop
    #créer la population initial à partir d'un dictionnaire
    def create_Dictionnairefinal(): 
            l1=list()
            l2=list()
            d=dict()
            with open('HSFinalSIVF.csv',  mode='r') as infile:
                    reader = csv.reader(infile, delimiter=';')
                    next(reader)
                    for rows in reader:
                        g = (rows[0],rows[1])                
                        l1.append(g) 
                        IS=(float(rows[2])/1000,float(rows[3]))  
                        l2.append(IS) 
                        
            DictF={k:v for k, v in zip(l1,l2)} # ==> size de dictF 891.558 
                                               
            return (DictF)
#----------------------------------------------------------------------------------------------------------------------------#
    def simInt_fsi(individu, dictgene): # moy_sim, moy_int et fitness par individu
        st = 0 # somme_taille = nbre d'interactions et de similarités par individu
        si = 0.0
        ss = 0.0
        alpha = 0.5 
        beta = 0.5 
        for i, item in enumerate(individu):
            #print("---------------------------")
            for j in range(i+1, len(individu)):
                g1 = individu[i]
                g2 = individu[j]
                #print("(g1, g2)=", g1,",",g2) 
                genes = (g1, g2)
                geness = (g2, g1)
                t = len(individu)
                v1 = dictgene.get(genes, 0) 
                v2 = dictgene.get(geness, 0) 
                
                if v1==0: 
                    v=v2 
                else: 
                    v=v1 
                if (v!=0):
                    
                    #print("       v=", v) # verif : RNASEH2A, ITPA ligne 862 969 cellule B862969                        
                    intr = v[0]
                    s = v[1]
                    si = si + intr 
                    ss = ss + s  
                    #st+=1
                    
                else: # calcul de GS2
                    intr = 0.0
                    s = 0.0 # GS2(g1, g2)
                    
            st = st + (t-i)
        #print("stt =", st)
        sint = si/st
        smoyen = ss/st
        
        f = alpha*(smoyen) + beta*(sint) # ---------F(S) = (W1*AVG_SIM) + (W2*AVG_Interaction)
        
        return sint, smoyen, f
        
    def fitness(candidates,args):
        print ('------------------------')
        print ("fitness")
        l_fit = []
        c = 0 
        for individu in candidates:
            
            sint, smoyen, f = simInt_fsi(individu, dictgene)
            #print("Int_moy =", sint)        
            #print("Sim_moy =", smoyen)
            print ('-- fitness =', f)
            l_fit.append(f)
            #c = c+1 # compteur pour itérer sur les individus       
            #print("community n°:",c,", size = ", len(individu))
            
        print("liste fitness =", l_fit)
        return (l_fit)     
#----------------------------------------------------------------------------------------------------------------------------#
    def population(candidates):
        l_map = []
        for ii in candidates:  # ii = liste de gènes = communauté
                l_map.append((len(ii),ii)) # ajoute (taille, communauté) à la liste l_map 
        #print(l_map)
        return l_map
#----------------------------------------------------------------------------------------------------------------------------#

    def Getbestgene(Score): # Score est une liste de couples (gène, score)
        listt = list()
        for i in Score:
            listt.append(i)
        bestgene = max(listt, key = itemgetter(1))[0] 
        scoreBestgene = max(listt, key = itemgetter(1))[1] #--donne moi le score de best_gene
        return (bestgene, scoreBestgene)  

    def GetWorstGene(Score): 
        listt = list()
        for i in Score:
            listt.append(i)
        Worstgene = min(listt, key = itemgetter(1))[0] 
        scoreWorstgene = min(listt, key = itemgetter(1))[1] 
        return (Worstgene , scoreWorstgene)  
    def get_som_si(index, candidate):
        si = 0.0
        ss = 0.0
        g1 = ''
        
        for j, item in enumerate(candidate):   
            g1=candidate[index] 
            g2=candidate[j]
            #print("(g1, g2)=", g1,",",g2)
            genes=(g1,g2) 
            geness=(g2,g1) 
            v1=dictgene.get(genes, 0) 
            v2=dictgene.get(geness,0) 
                        
            if v1==0: 
                v=v2 
            else:
                v=v1 
                                
            if (v!=0):
                #print("       v=", v)
                i = float(v[0])    
                s = float(v[1]) 
            else: 
                i=0.0
                s=0.0 # s = GS2(g1, g2)
            si=si+i
            ss=ss+s
            
        return g1, si, ss

    def Getscorepergene(candidate):    
            g1=''
            lscore=list()
            for index , cs in enumerate(candidate): 
                #print("---------------------------")
                si = 0.0
                ss = 0.0
                t = 0
                g1, si, ss = get_som_si(index, candidate)
                t = len(candidate)
                
                #print("t =", t)
                #avg_sim = ss / (2*(t-1))
                #avg_int = si / (2*(t-1))
                
                avg_sim = ss / (t-1)
                avg_int = si / (t-1)
                score = avg_sim + avg_int
                lscore.append((g1, score))
            print ("score par gene :", lscore)
            return(lscore)


    def GetGeneInteraction(mutant, bestgene):
        l = list()
        for k in keyssample:
            if bestgene == k[0]:
                l.append(k[1])
            if bestgene == k[1]:
                l.append(k[0])
        l = list(set(l))       
        for i, item in enumerate(l):
            if(l[i] in mutant):
                #print("l[i] =", l[i])
                l.remove(l[i])
        #print("taille liste après supp:", len(l)) # l=liste des gènes distincts qui interagissent avec bestgene et ne figurent pas dans mutant
        if (not l): # score_best = 0.0
            choix = rand.choice([gene for gene in mutant if gene != bestgene])
        else:
            choix = rand.choice(l)
            
        print("gene interact :", choix)
        return choix
#----------------------------------------------------------------------------------------------------------------------------#

    def mutate_gene (random, child, args):
            print ("--------------mutation OCM1 -------------")
            m=list()
            bounder = args['_ec'].bounder
            seuil=0.5
            m=copy.copy(child) # signification -- après croisement, on copie les 2 enfants pour éviter tt changement
            mutant= list(set(m)) # set(m) : No duplicate members => Pour éliminer la redondance dans m 
            
            # crossover offsprings dans mutation = S : Individu Enfant obtenu après le croisement
            print (" crossover offsprings dans mutation = ", mutant , type(mutant)) 
            print ("Size mutant =", len(mutant))
            s1 = Getscorepergene(mutant)
            worst, sw = GetWorstGene(s1) # sw = score_worst
            best,sb = Getbestgene(s1) # sb = score_best
            print ("Best =", best)
            print ("Score_Best =", sb)
            print ("Worst =", worst)
            print ("Score_Worst =", sw)
            GeneIn = GetGeneInteraction(mutant, best)
            
            for i, m in enumerate(mutant): # m = gene dans mutant
                 

                    gene = mutant[i].strip("' '")
                    if GeneIn=='None':
                        pass
                    if gene != GeneIn:
                        if  (seuil <= sw): 
                            #print 'gene= ',gene
                            print ("***** Insertion *****")
                            mutant.insert(len(mutant)+1, GeneIn)
                            mutant = list(mutant) # pour MAJ la taille de mutant
                            candidate = copy.copy(mutant)
                            print ("Offspring after mutate =", candidate) # S' : Individu enfant issu du processus de mutation.
                            print ("Fin mutation ")
                            logger.debug("Finished with child num after mutate {0} ".format(candidate) )
                            # This message should go to the log file

                            return candidate

                        else :
                            for k,p in enumerate(mutant):
                                if (mutant[k]) == worst :
                                    print ("Replace")
                                    mutant[k] = str(GeneIn)
                                    mutant=list(mutant) # MAJ la taille
                                    candidate = copy.copy(mutant)
                                    #candidate.remove('None')
                                    print ("offspring after mutate ",candidate)
                                    print ("Fin mutation ")
                                    logger.debug("Finished with child num after mutate {0} ".format(candidate) )
                                    return candidate    

#-------------------------------------------------------------------------------------------------------
#generation de fichier excl pour sauvegarder les données

    def getInfos(taille, best_fit, duration, best_candidate, dictgene):
    
        sint, smoy, f = simInt_fsi(best_candidate, dictgene)
        
        
        df = pd.DataFrame([[taille, '', best_fit, duration, sint, smoy, str(best_candidate)]], columns=['Taille', '% recouvrement', 'f_si','Duration', 'Int_moy', 'Sim_moy', 'Individu'])
        return df

    def generate_file(file, df):
        
        if os.path.isfile(file):  # if file already exists append to existing file
            workbook = openpyxl.load_workbook(file)  # load workbook if already exists
            sheet = workbook['my_sheet_si']  # declare the active sheet 
    
            # append the dataframe results to the current excel file
            for row in dataframe_to_rows(df, header = False, index = False):
                sheet.append(row)
            workbook.save(file)  # save workbook
            workbook.close()  # close workbook
        else:  # create the excel file if doesn't already exist
            with pd.ExcelWriter(path = file, engine = 'openpyxl') as writer:
                df.to_excel(writer, index = False, sheet_name = 'my_sheet_si')
        print("data appended.")

#-------------------------------------------------------------------------------------------

    def testObservers(population, num_generations, num_evaluations, args):
        f=open('bestSamplealphai.txt', 'a')
        f1=open('samples.txt', 'a')
        logger.debug("------ ----- ----- Generation num ----- ----- ----- {0}".format(num_generations) )
        print ('------ ----- ----- Generation num ------ ----- ----', num_generations)
        print ("----------- population num ", num_generations)

        for p in population:

                logger.debug(p)   
        inspyred.ec.observers.stats_observer(population, num_generations, num_evaluations, args) 
        #Print the statistics of the evolutionary computation to the screen.
        inspyred.ec.observers.plot_observer(population, num_generations, num_evaluations, args)          

    #dfSI=create_DataFrame() 
    dictgene = create_Dictionnairefinal()
    keys = dictgene.keys()
    values = dictgene.values()
    keyssample = dictgene.keys() # les cles de dictio
    rand = Random()
    
    algorithm = inspyred.ec.GA(rand)
    algorithm.terminator = inspyred.ec.terminators.generation_termination
    algorithm.observer = testObservers 
    #algorithm.observer = inspyred.ec.observers.best_observer
    algorithm.selector = inspyred.ec.selectors.tournament_selection
    algorithm.replacer = inspyred.ec.replacers. steady_state_replacement 
    '''#steady_state_replacement les chromosomes les plus mauvais sont retires et remplacer par les nouveaux. 
    En consequence, le reste de la 		population survit a la nouvelle generation.  '''
    
    algorithm.variator = [inspyred.ec.variators.n_point_crossover, inspyred.ec.variators.mutator(mutate_gene)]
    algorithm.archiver = inspyred.ec.archivers.population_archiver
    algorithm.migrator  = inspyred.ec.migrators.default_migration
    #inspyred.ec.analysis.generation_plot(file_handler, errorbars=True)
    
    final_pop = algorithm.evolve(generator = generate_pop,
                         evaluator = fitness,
                         pop_size = 30, 
                         maximize = True, # maximise la fct fitness
                         num_selected = 2, 
                         tournament_size = 10, 
                         max_evaluations = 1000, 
                         mutation_rate = 0.01, # 0.1 -> 0.5
                         #num_crossover_points=2, #pour appliquer le croisement 2points
                         crossover_rate = 0.8, #  0.5 -> 0.8 
                         max_generations = 2 # 50 -> 100
                         )
        
    if display:
            
            fp = list(final_pop)
            best = max(final_pop)    
            print ('Pop final: \n{0}'.format(str(fp)))
            print('Best Solution: \n{0}'.format(str(best)))
            logger.debug('Population final: \n{0}'.format(str(fp))) #le prob ici il fait debug plusieurs fois la pop
            logger.debug('Best Solution: {0}'.format(str(best)))
            end_time = datetime.now()
            print('Duration: {}'.format(end_time - start_time))
            
            duration = str(end_time - start_time)
            taille = len(best.candidate)
            
            df = getInfos(taille, best.fitness, duration, best.candidate, dictgene)#créer un fichier excl au fur et à mesure 
            file = 'test_si.xlsx'
            generate_file(file, df)
            
if __name__ == '__main__':
    print ("debut \n")
    print("Pop Initial \n")
    main(display=True)