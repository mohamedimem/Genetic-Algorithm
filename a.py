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

def setup_logger(metric_name):
    """Setup logger for each metric type"""
    logger = logging.getLogger(f'inspyred.ec.{metric_name}')
    logger.propagate = False 
    logger.setLevel(logging.DEBUG)
    file_handler = logging.handlers.RotatingFileHandler(f'inspyred_{metric_name}.log', mode='w')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger

def main_algorithm(location_metric='jaccard_max', run_number=1, display=False):
    """
    Main genetic algorithm function
    location_metric: 'jaccard_max', 'adamic_max', or 'gaussian_max'
    run_number: current run number (1-20)
    """
    
    start_time = datetime.now()
    rand = Random() 
    tabmax = 40

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
            print(f"Run {run_number} - {location_metric}: {pop}")
        return pop
    
    def create_Dictionnairefinal(): 
        l1=list()
        l2=list()
        d=dict()
        
        # Define column mapping for different metrics
        metric_column = {
            'jaccard_max': 4,
            'adamic_max': 5, 
            'gaussian_max': 6
        }
        
        with open('HSFinalSIVF.csv', mode='r') as infile:
            reader = csv.reader(infile, delimiter=';')
            next(reader)
            for rows in reader:
                g = (rows[0],rows[1])                
                l1.append(g) 
                
                interaction = float(rows[2])/1000
                similarity = float(rows[3])
                
                # Handle the specified location metric
                col_index = metric_column[location_metric]
                if len(rows) > col_index and rows[col_index] and rows[col_index].strip():
                    location_value = float(rows[col_index])
                    has_location_data = True
                else:
                    location_value = None
                    has_location_data = False
                
                IS = (interaction, similarity, location_value, has_location_data)
                l2.append(IS) 
                
        DictF = {k:v for k, v in zip(l1,l2)}
        return DictF

    def simInt_fsi(individu, dictgene):
        st = 0
        si = 0.0
        ss = 0.0
        sl = 0.0
        location_count = 0
        has_any_location = False
        
        alpha = 0.4
        beta = 0.4
        gamma = 0.2
        
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
                    location_val = v[2]
                    has_data = v[3]
                    
                    si = si + intr 
                    ss = ss + s  
                    
                    if has_data:
                        sl = sl + location_val
                        location_count += 1
                        if location_val > 0:
                            has_any_location = True
                    
                else:
                    intr = 0.0
                    s = 0.0
                    
            st = st + (t-i)
        
        sint = si/st
        smoyen = ss/st
        
        if location_count > 0 and has_any_location:
            slocation = sl/location_count
            f = alpha*(smoyen) + beta*(sint) + gamma*(slocation)
        else:
            f = 0.5*(smoyen) + 0.5*(sint)
        
        return sint, smoyen, f
        
    def fitness(candidates,args):
        if display:
            print(f'Run {run_number} - {location_metric} - fitness calculation')
        l_fit = []
        for individu in candidates:
            sint, smoyen, f = simInt_fsi(individu, dictgene)
            if display:
                print(f'-- fitness = {f}')
            l_fit.append(f)
        return l_fit

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
        sl = 0.0
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
                if v[3]:  # has location data
                    l = float(v[2])
                else:
                    l = 0.0
            else: 
                i=0.0
                s=0.0
                l=0.0
            si=si+i
            ss=ss+s
            sl=sl+l
            
        return g1, si, ss, sl

    def Getscorepergene(candidate):    
        g1=''
        lscore=list()
        for index , cs in enumerate(candidate): 
            si = 0.0
            ss = 0.0
            sl = 0.0
            t = 0
            g1, si, ss, sl = get_som_si(index, candidate)
            t = len(candidate)
            
            avg_sim = ss / (t-1)
            avg_int = si / (t-1)
            avg_loc = sl / (t-1)
            score = avg_sim + avg_int + avg_loc
            lscore.append((g1, score))
        if display:
            print(f"Run {run_number} - {location_metric} - score par gene:", lscore)
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
            print(f"Run {run_number} - {location_metric} - gene interact:", choix)
        return choix

    def mutate_gene(random, child, args):
        if display:
            print(f"Run {run_number} - {location_metric} - mutation OCM1")
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
                if (seuil <= sw): 
                    mutant.insert(len(mutant)+1, GeneIn)
                    mutant = list(mutant)
                    candidate = copy.copy(mutant)
                    return candidate
                else :
                    for k,p in enumerate(mutant):
                        if (mutant[k]) == worst :
                            mutant[k] = str(GeneIn)
                            mutant=list(mutant)
                            candidate = copy.copy(mutant)
                            return candidate    

    def log_result_to_file(result, results_file):
        """Log result to a simple text file"""
        with open(results_file, 'a') as f:
            f.write(f"Run {result['run']} - {result['metric'].upper()}\n")
            f.write(f"Fitness: {result['fitness']:.6f}\n")
            f.write(f"Duration: {result['duration']}\n")
            f.write(f"Community Size: {result['size']}\n")
            f.write(f"Best Solution: {', '.join(result['candidate'])}\n")
            f.write("-" * 80 + "\n\n")

    def testObservers(population, num_generations, num_evaluations, args):
        if display and num_generations % 10 == 0:  # Print every 10 generations
            print(f'Run {run_number} - {location_metric} - Generation {num_generations}')
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
    algorithm.migrator = inspyred.ec.migrators.default_migration
    
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
        print(f'Run {run_number} - {location_metric} - Best Solution: {best.candidate} : {best.fitness}')
        print(f'Duration: {duration}')
    
    # Return results for logging
    return {
        'run': run_number,
        'metric': location_metric,
        'fitness': best.fitness,
        'duration': duration,
        'candidate': best.candidate,
        'size': taille
    }

def log_result_to_file(result, results_file):
    """Log result to a simple text file"""
    with open(results_file, 'a') as f:
        f.write(f"Run {result['run']} - {result['metric'].upper()}\n")
        f.write(f"Fitness: {result['fitness']:.6f}\n")
        f.write(f"Duration: {result['duration']}\n")
        f.write(f"Community Size: {result['size']}\n")
        f.write(f"Best Solution: {', '.join(result['candidate'])}\n")
        f.write("-" * 80 + "\n\n")

def run_experiments():
    """Run all experiments: 20 runs x 3 metrics = 60 total runs"""
    
    print("=== Starting Multi-Run Genetic Algorithm Experiments ===")
    print("3 metrics × 20 runs each = 60 total runs")
    print()
    
    metrics = ['jaccard_max', 'adamic_max', 'gaussian_max']
    all_results = []
    
    # Create results file (simple text file)
    results_file = f'GA_MultiRun_Results_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt'
    
    # Write header to results file
    with open(results_file, 'w') as f:
        f.write("GENETIC ALGORITHM MULTI-RUN EXPERIMENT RESULTS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Experiments: 60 (20 runs × 3 metrics)\n")
        f.write("=" * 80 + "\n\n")
    
    for metric in metrics:
        print(f"\n--- Starting 20 runs for {metric.upper()} ---")
        
        # Setup logger for this metric
        logger = setup_logger(metric)
        
        metric_results = []
        
        # Write metric section header
        with open(results_file, 'a') as f:
            f.write(f"\n{'='*20} {metric.upper()} RESULTS {'='*20}\n\n")
        
        for run in range(1, 21):  # 20 runs
            print(f"Running {metric} - Run {run}/20...")
            
            try:
                result = main_algorithm(
                    location_metric=metric, 
                    run_number=run, 
                    display=(run == 1)  # Only show details for first run of each metric
                )
                
                metric_results.append(result)
                all_results.append(result)
                
                # Log to individual metric log file
                logger.info(f"Run {run} - Fitness: {result['fitness']:.6f} - Size: {result['size']} - Duration: {result['duration']}")
                
                # Log to main results file
                log_result_to_file(result, results_file)
                
                print(f"  ✓ Run {run} completed - Fitness: {result['fitness']:.6f}")
                
            except Exception as e:
                print(f"  ✗ Error in {metric} run {run}: {e}")
                logger.error(f"Run {run} failed: {e}")
                
                # Log error to results file
                with open(results_file, 'a') as f:
                    f.write(f"Run {run} - {metric.upper()} - ERROR: {e}\n")
                    f.write("-" * 80 + "\n\n")
        
        # Summary for this metric
        if metric_results:
            best_run = max(metric_results, key=lambda x: x['fitness'])
            worst_run = min(metric_results, key=lambda x: x['fitness'])
            avg_fitness = sum(r['fitness'] for r in metric_results) / len(metric_results)
            
            summary_text = f"""
{metric.upper()} SUMMARY:
  Completed runs: {len(metric_results)}/20
  Best fitness:   {best_run['fitness']:.6f} (Run {best_run['run']})
  Worst fitness:  {worst_run['fitness']:.6f} (Run {worst_run['run']})
  Average fitness: {avg_fitness:.6f}
  Best solution size: {best_run['size']} genes
  Best solution: {', '.join(best_run['candidate'])}
"""
            
            print(summary_text)
            
            # Write summary to results file
            with open(results_file, 'a') as f:
                f.write(summary_text + "\n")
                f.write("=" * 80 + "\n\n")
            
            # Log summary to individual log
            logger.info(f"SUMMARY - Best: {best_run['fitness']:.6f}, Worst: {worst_run['fitness']:.6f}, Avg: {avg_fitness:.6f}")
    
    # Overall summary
    summary_text = f"""
=== OVERALL EXPERIMENT SUMMARY ===
Total runs completed: {len(all_results)}/60
Completion time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Results saved to: {results_file}
"""
    
    print(summary_text)
    
    if all_results:
        overall_best = max(all_results, key=lambda x: x['fitness'])
        
        final_summary = f"""
BEST OVERALL SOLUTION:
  Metric: {overall_best['metric'].upper()}
  Run: {overall_best['run']}
  Fitness: {overall_best['fitness']:.6f}
  Community Size: {overall_best['size']} genes
  Duration: {overall_best['duration']}
  Solution: {', '.join(overall_best['candidate'])}
"""
        
        print(final_summary)
        
        # Write final summary to results file
        with open(results_file, 'a') as f:
            f.write(summary_text)
            f.write(final_summary)
            f.write("\n" + "=" * 80 + "\n")
            f.write("EXPERIMENT COMPLETED\n")
    
    return all_results, results_file

if __name__ == '__main__':
    results, filename = run_experiments()
    print(f"\nExperiment completed! Check {filename} for detailed results.")