import matplotlib.pyplot as plt
from numpy import*
import numpy as np
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
from operator import itemgetter
import itertools
import random
from datetime import datetime
from multiprocessing import Manager 

import logging
import logging.handlers

def setup_logger(run_number):
    """Setup logger for each run"""
    logger = logging.getLogger(f'inspyred.ec.run_{run_number}')
    logger.propagate = False 
    logger.setLevel(logging.DEBUG)
    file_handler = logging.handlers.RotatingFileHandler(f'inspyred_run_{run_number}.log', mode='w')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger

def main_single_run(run_number=1, display=False):
    """Single run of the original algorithm"""
    
    start_time = datetime.now()
    rand = Random() 
    tabmax = 40

    # Setup logger for this run
    logger = setup_logger(run_number)

    def getsizeCommunitie(): 
        var = rand.randint(5,tabmax)
        return var
    
    def getGeneFromFile():
        f2 = open ('pathW2.txt','r')
        l = list() 
        
        for  i in  f2.readlines(): 
            lignes = i.strip('\n ') 
            l.append(lignes)
        
        nv = list()
        v = getsizeCommunitie()
        
        for i in range(0,v):
            i = rand.choice(l) 
            l.remove(i) 
            nv.append(i)
        return nv       
    
    def generate_pop(random, args): 
        pop = getGeneFromFile()
        if display:
            print(f"Run {run_number}: {pop}")
        return pop
    
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
                
        DictF={k:v for k, v in zip(l1,l2)}
        return (DictF)

    def simInt_fsi(individu, dictgene):
        st = 0
        si = 0.0
        ss = 0.0
        alpha = 0.5 
        beta = 0.5 
        for i, item in enumerate(individu):
            for j in range(i+1, len(individu)):
                g1 = individu[i]
                g2 = individu[j]
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
                    intr = v[0]
                    s = v[1]
                    si = si + intr 
                    ss = ss + s  
                else:
                    intr = 0.0
                    s = 0.0
                    
            st = st + (t-i)
        
        sint = si/st
        smoyen = ss/st
        f = alpha*(smoyen) + beta*(sint)
        
        return sint, smoyen, f
        
    def fitness(candidates,args):
        if display:
            print(f'Run {run_number} - fitness calculation')
        l_fit = []
        for individu in candidates:
            sint, smoyen, f = simInt_fsi(individu, dictgene)
            if display:
                print(f'-- fitness = {f}')
            l_fit.append(f)
        return (l_fit)

    def Getbestgene(Score):
        listt = list()
        for i in Score:
            listt.append(i)
        bestgene = max(listt, key = itemgetter(1))[0] 
        scoreBestgene = max(listt, key = itemgetter(1))[1]
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
            genes=(g1,g2) 
            geness=(g2,g1) 
            v1=dictgene.get(genes, 0) 
            v2=dictgene.get(geness,0) 
                        
            if v1==0: 
                v=v2 
            else:
                v=v1 
                                
            if (v!=0):
                i = float(v[0])    
                s = float(v[1]) 
            else: 
                i=0.0
                s=0.0
            si=si+i
            ss=ss+s
            
        return g1, si, ss

    def Getscorepergene(candidate):    
        g1=''
        lscore=list()
        for index , cs in enumerate(candidate): 
            si = 0.0
            ss = 0.0
            t = 0
            g1, si, ss = get_som_si(index, candidate)
            t = len(candidate)
            
            avg_sim = ss / (t-1)
            avg_int = si / (t-1)
            score = avg_sim + avg_int
            lscore.append((g1, score))
        if display:
            print(f"Run {run_number} - score par gene:", lscore)
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
                l.remove(l[i])
        if (not l):
            choix = rand.choice([gene for gene in mutant if gene != bestgene])
        else:
            choix = rand.choice(l)
            
        if display:
            print(f"Run {run_number} - gene interact:", choix)
        return choix

    def mutate_gene (random, child, args):
        if display:
            print(f"Run {run_number} - mutation OCM1")
        m=list()
        bounder = args['_ec'].bounder
        seuil=0.5
        m=copy.copy(child)
        mutant= list(set(m))
        
        s1 = Getscorepergene(mutant)
        worst, sw = GetWorstGene(s1)
        best,sb = Getbestgene(s1)
        GeneIn = GetGeneInteraction(mutant, best)
        
        for i, m in enumerate(mutant):
            gene = mutant[i].strip("' '")
            if GeneIn=='None':
                pass
            if gene != GeneIn:
                if  (seuil <= sw): 
                    mutant.insert(len(mutant)+1, GeneIn)
                    mutant = list(mutant)
                    candidate = copy.copy(mutant)
                    logger.debug("Finished with child num after mutate {0} ".format(candidate) )
                    return candidate
                else :
                    for k,p in enumerate(mutant):
                        if (mutant[k]) == worst :
                            mutant[k] = str(GeneIn)
                            mutant=list(mutant)
                            candidate = copy.copy(mutant)
                            logger.debug("Finished with child num after mutate {0} ".format(candidate) )
                            return candidate    

    def testObservers(population, num_generations, num_evaluations, args):
        logger.debug("------ ----- ----- Generation num ----- ----- ----- {0}".format(num_generations) )
        if display and num_generations % 10 == 0:  # Print every 10 generations
            print(f'Run {run_number} - Generation {num_generations}')
            for p in population:
                logger.debug(p)   
            inspyred.ec.observers.stats_observer(population, num_generations, num_evaluations, args)

    # Main algorithm execution
    dictgene = create_Dictionnairefinal()
    keys = dictgene.keys()
    values = dictgene.values()
    keyssample = dictgene.keys()
    rand = Random()
    
    algorithm = inspyred.ec.GA(rand)
    algorithm.terminator = inspyred.ec.terminators.generation_termination
    algorithm.observer = testObservers 
    algorithm.selector = inspyred.ec.selectors.tournament_selection
    algorithm.replacer = inspyred.ec.replacers.steady_state_replacement 
    algorithm.variator = [inspyred.ec.variators.n_point_crossover, inspyred.ec.variators.mutator(mutate_gene)]
    algorithm.archiver = inspyred.ec.archivers.population_archiver
    algorithm.migrator  = inspyred.ec.migrators.default_migration
    
    final_pop = algorithm.evolve(generator = generate_pop,
                         evaluator = fitness,
                         pop_size = 30, 
                         maximize = True,
                         num_selected = 2, 
                         tournament_size = 10, 
                         max_evaluations = 1000, 
                         mutation_rate = 0.01,
                         crossover_rate = 0.8,
                         max_generations = 100
                         )
        
    best = max(final_pop)    
    end_time = datetime.now()
    duration = str(end_time - start_time)
    taille = len(best.candidate)
    
    if display:
        print(f'Run {run_number} - Best Solution: {best.candidate} : {best.fitness}')
        print(f'Duration: {duration}')
    
    # Log to main logger
    logger.info(f"Run {run_number} COMPLETED - Fitness: {best.fitness:.6f} - Size: {taille} - Duration: {duration}")
    
    # Return results
    return {
        'run': run_number,
        'fitness': best.fitness,
        'duration': duration,
        'candidate': best.candidate,
        'size': taille,
        'avg_similarity': simInt_fsi(best.candidate, dictgene)[1],
        'avg_interaction': simInt_fsi(best.candidate, dictgene)[0]
    }

def log_result_to_file(result, results_file):
    """Log result to a simple text file in tabular format"""
    with open(results_file, 'a') as f:
        f.write(f"{result['run']:2d} | {result['fitness']:8.6f} | {result['size']:2d} | {result['avg_similarity']:8.6f} | {result['avg_interaction']:8.6f} | {result['duration']:>12} | {', '.join(result['candidate'])}\n")

def run_20_experiments():
    """Run original algorithm 20 times"""
    
    print("=== Starting 20 Runs of Original Algorithm ===")
    print("Using original fitness function: 0.5*similarity + 0.5*interaction")
    print()
    
    all_results = []
    
    # Create results file
    results_file = f'Original_Algorithm_20Runs_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt'
    
    # Write header
    with open(results_file, 'w') as f:
        f.write("ORIGINAL GENETIC ALGORITHM - 20 RUNS EXPERIMENT RESULTS\n")
        f.write("=" * 120 + "\n")
        f.write(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("Fitness Function: 0.5 * AVG_SIMILARITY + 0.5 * AVG_INTERACTION\n")
        f.write("Population Size: 30 | Max Generations: 100 | Community Size: 5-40 genes (random)\n")
        f.write("=" * 120 + "\n\n")
        f.write("RUN RESULTS TABLE:\n")
        f.write("-" * 120 + "\n")
        f.write("RUN | FITNESS  | SIZE | AVG_SIM  | AVG_INT  |   DURATION   | SOLUTION (GENES)\n")
        f.write("-" * 120 + "\n")
    
    for run in range(1, 21):  # 20 runs
        print(f"Running Original Algorithm - Run {run}/20...")
        
        try:
            result = main_single_run(
                run_number=run, 
                display=(run == 1)  # Only show details for first run
            )
            
            all_results.append(result)
            
            # Log to results file in table format
            log_result_to_file(result, results_file)
            
            print(f"  ✓ Run {run} completed - Fitness: {result['fitness']:.6f} - Size: {result['size']} genes")
            
        except Exception as e:
            print(f"  ✗ Error in run {run}: {e}")
            
            # Log error to results file
            with open(results_file, 'a') as f:
                f.write(f"{run:2d} | ERROR: {str(e)[:80]}\n")
    
    # Calculate summary statistics
    if all_results:
        best_run = max(all_results, key=lambda x: x['fitness'])
        worst_run = min(all_results, key=lambda x: x['fitness'])
        avg_fitness = sum(r['fitness'] for r in all_results) / len(all_results)
        avg_size = sum(r['size'] for r in all_results) / len(all_results)
        
        # Sort by fitness for analysis
        sorted_results = sorted(all_results, key=lambda x: x['fitness'], reverse=True)
        
        summary_text = f"""
{"-" * 120}

ANALYSIS SUMMARY:
{"-" * 40}
Completed runs: {len(all_results)}/20
Best fitness:   {best_run['fitness']:.6f} (Run {best_run['run']}, Size: {best_run['size']} genes)
Worst fitness:  {worst_run['fitness']:.6f} (Run {worst_run['run']}, Size: {worst_run['size']} genes)
Average fitness: {avg_fitness:.6f}
Average community size: {avg_size:.1f} genes

FITNESS RANKING (Top 10):
{"-" * 40}
"""
        
        # Add top 10 results for analysis
        for i, result in enumerate(sorted_results[:10], 1):
            summary_text += f"{i:2d}. Run {result['run']:2d} - Fitness: {result['fitness']:.6f} - Size: {result['size']:2d} genes\n"
        
        summary_text += f"\nSIZE DISTRIBUTION ANALYSIS:\n{'-' * 40}\n"
        
        # Size analysis
        size_counts = {}
        for result in all_results:
            size = result['size']
            if size not in size_counts:
                size_counts[size] = []
            size_counts[size].append(result['fitness'])
        
        # Sort by size
        for size in sorted(size_counts.keys()):
            fitness_list = size_counts[size]
            avg_fitness_for_size = sum(fitness_list) / len(fitness_list)
            max_fitness_for_size = max(fitness_list)
            summary_text += f"Size {size:2d}: {len(fitness_list)} runs, Avg Fitness: {avg_fitness_for_size:.6f}, Best: {max_fitness_for_size:.6f}\n"
        
        summary_text += f"\nCOMPONENT ANALYSIS:\n{'-' * 40}\n"
        summary_text += f"Best Run ({best_run['run']}) Details:\n"
        summary_text += f"  Similarity Component: {best_run['avg_similarity']:.6f} (weight: 0.5)\n"
        summary_text += f"  Interaction Component: {best_run['avg_interaction']:.6f} (weight: 0.5)\n"
        summary_text += f"  Combined Fitness: {best_run['fitness']:.6f}\n"
        summary_text += f"  Genes: {', '.join(best_run['candidate'])}\n"
        
        print(summary_text)
        
        # Write summary to results file
        with open(results_file, 'a') as f:
            f.write(summary_text + "\n")
            f.write("=" * 120 + "\n")
            f.write(f"Experiment completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    print(f"\nExperiment completed! Results saved to: {results_file}")
    return all_results, results_file

if __name__ == '__main__':
    print("Starting Original Algorithm - 20 Runs Experiment\n")
    results, filename = run_20_experiments()
    print(f"\nAll runs completed! Check {filename} for detailed results.")