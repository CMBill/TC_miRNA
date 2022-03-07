# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:02:34 2021

@author: 12870
"""

import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
import copy
import random
import numpy as np

Init_score = -1
data= pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\generank_selected.csv')
sample_name = data.columns.tolist()


class Life(object):
    def __init__(self, aGene=None):
        self.gene = aGene
        self.score = Init_score


class GA(object):

    def __init__(self,
                 CrossRate,
                 MutationRage,
                 individual,
                 genenumber,
                 aMatchFun=lambda life: 1):
        self.crossRate = CrossRate
        self.mutationRate = MutationRage
        self.popsize = individual
        self.chromlength = genenumber
        self.matchFun = aMatchFun
        self.pop = []
        self.best = Life(np.random.randint(0, 1, self.chromlength))

        self.gene = np.random.randint(0, 1, self.chromlength)
        self.score = -1

        self.generation = 0
        self.crossCount = 0
        self.mutationCount = 0
        self.bounds = 0.0
        self.geneEncoding()

    def geneEncoding(self):
        self.pop = []
        count = 0
        while count < self.popsize:
            temp = []
            has = False
            for j in range(self.chromlength):
                gene = random.randint(0, 1)
                if gene == 1:
                    has = True
                temp.append(gene)
            life = Life(temp)
            random.shuffle(temp)
            if has:
                self.pop.append(life)
                count += 1

    def calfitness(self):
        self.bounds = 0.0
        self.best.score = copy.deepcopy(self.score)
        self.best.gene = copy.deepcopy(self.gene)
        for life in self.pop:
            life.score = self.matchFun(life)
            self.bounds += life.score
            if self.best.score < life.score:
                self.best = life

        if self.score < self.best.score:
            self.score = copy.deepcopy(self.best.score)
            self.gene = copy.deepcopy(self.best.gene)


    def crossover(self, parent1, parent2):
        cpoint1 = random.randint(0, self.chromlength - 1)  
        cpoint2 = random.randint(cpoint1, self.chromlength - 1)  

        for point in range(len(parent1.gene)):
            if (point >= cpoint1) and (point <= cpoint2):
                parent1.gene[point], parent2.gene[point] = parent2.gene[point], parent1.gene[point]

        self.crossCount += 1
        return parent1.gene

    def mutation(self, gene):
       
        mpoint = random.randint(0, self.chromlength - 1)
        newGene = gene[:]  
        if newGene[mpoint] == 1:
               newGene[mpoint] = 0
        else:
            newGene[mpoint]= 1
        self.mutationCount += 1
        return newGene

    def getOne(self):
        
        r = random.uniform(0, self.bounds)
        for life in self.pop:
            r -= life.score
            if r <= 0:
                return life

        raise Exception("wrong", self.bounds)

    def newChild(self):
        
        parent1 = self.getOne()
        rate1 = random.random()

        if rate1 < self.crossRate:
           
            parent2 = self.getOne()
            gene = self.crossover(parent1, parent2)
        else:
            gene = parent1.gene

        rate2 = random.random()
        if rate2 < self.mutationRate:
            gene = self.mutation(gene)

        return Life(gene)

    def next(self):
        self.calfitness()
        newpop = []
        newpop.append(self.best)  
        newpop[0].gene = copy.deepcopy(self.gene)
        newpop[0].score = copy.deepcopy(self.score)
        while len(newpop) < self.popsize:
            newpop.append(self.newChild())
        self.pop = newpop
        self.generation += 1


class FeatureSelection(object):
    def __init__(self, individual=20):
        
        self.columns = sample_name
        self.train_data = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.xlsx')
                       
        self.popsize = individual
        self.ga = GA(CrossRate=0.7,
                     MutationRage=0.1,
                     individual=self.popsize,
                     genenumber=len(self.columns) - 1,
                     aMatchFun=self.matchFun())
        
    def acc_score(self, gene):
        print(gene)
        features = self.columns[1:]
        features_name = []
        for i in range(len(gene)):
            if gene[i] == 1:
                features_name.append(features[i])

        labels = np.array(self.train_data['Type'], dtype=int)      
        clf = LogisticRegression(penalty='l2', random_state=0, C=1)
        score = cross_val_score(clf,self.train_data[features_name], labels, cv=10).mean()  # 10次交叉验证
        return score
    
    def matchFun(self):
        return lambda life: self.acc_score(life.gene)
    
    def run(self, n=0):
        acc_list = []
        iteration = [index for index in range(1, n + 1)]
        while n > 0:
            self.ga.next()
            acc = self.ga.score                     
            acc_list.append(acc)
            print(("第%d代 : 当前最好特征组合ACC：%f") % (self.ga.generation, acc))
            n -= 1  

        print('当前最好特征组合:')
        string = []
        flag = 0
        features = self.columns[1:]
        for index in self.ga.gene:                                 
            if index == 1:
                string.append(features[flag])
            flag += 1
        print(string)
        print('最高acc：', self.ga.score)                     

        ''''''
        plt.plot(iteration, acc_list)
        plt.xlabel('Generation')
        plt.ylabel('ACC(%)')
        plt.show()
        print('Each generation of optimality',acc_list)


def main():
    fs = FeatureSelection(individual=50)
    rounds =120    
    fs.run(rounds)


if __name__ == '__main__':
    main()
    

