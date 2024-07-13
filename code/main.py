from time import time

import numpy as np
import os
import pandas as pd
import warnings

from algorithm.BaseMissing import BaseMissing
from algorithm.SCIILP import SCIILP
from algorithm.SCILP import SCILP
from algorithm.SCILN import SCILN
from util.Assist import Assist
from util.DataHandler import DataHandler as dh
from util.FileHandler import FileHandler as fh

import logging

enable_debug = False
stop_check = False

logging.basicConfig(level=logging.INFO)
if enable_debug:
    logging.getLogger().setLevel(logging.DEBUG)
else:
    logging.getLogger().setLevel(logging.INFO)

warnings.filterwarnings("ignore")


class CompTest:
    def __init__(self, filename, label_idx=-1, exp_cluster=None, header=None, index_col=None):
        self.name = filename.split('/')[-1].split('.')[0]
        self.REPEATLEN = 5
        self.RBEGIN = 998
        self.exp_cluster = exp_cluster
        self.label_idx = label_idx
        self.fh = fh()
        self.db = self.fh.readCompData(filename, label_idx=label_idx, header=header, index_col=index_col)
        self.dh = dh()
        self.delCompRowIndexList = self.fh.findMis()
        self.totalsize = self.db.shape[0]
        self.size = self.totalsize * 1

        self.ALGNUM = 3
        self.totalTime = np.zeros((1, self.ALGNUM))
        self.purityAlg = np.zeros((self.ALGNUM))
        self.riAlg = np.zeros((self.ALGNUM))
        self.ariAlg = np.zeros((self.ALGNUM))
        self.fmeasureAlg = np.zeros((self.ALGNUM))
        self.totalpurity = np.zeros((1, self.ALGNUM))
        self.totalfmeasure = np.zeros((1, self.ALGNUM))
        self.totalri = np.zeros((1, self.ALGNUM))
        self.totalari = np.zeros((1, self.ALGNUM))

        self.alg_flags = [1, 1, 1]

        self.K_SCI = exp_cluster
        self.KDistance = 5
        self.K_Candidate_SCI = 1
        self.epsilon_SCI = 0.01

    def SCI_LN(self):
        print("SCI_LN begin!")
        algIndex = 0
        if self.alg_flags[algIndex]:
            startTime = time()
            for repeat in range(self.RBEGIN, self.RBEGIN + self.REPEATLEN):
                seed = 100 + self.size * 1011 + repeat * 811
                cells, mask = self.dh.getMisCelMul(self.db, self.label_idx, self.size)

                SCI = SCILN(self.db, self.label_idx, self.exp_cluster, cells, mask)
                SCI.setK(self.K_SCI)
                SCI.setK_Candidate(self.K_Candidate_SCI)
                SCI.setEpsilon(self.epsilon_SCI)
                SCI.setKDistance(self.KDistance)
                SCI.setDelCompRowIndexList([])
                SCI.mainSCI()

                cluster_dict = {}
                for i in range(len(SCI.clusterMembersList)):
                    for j in SCI.clusterMembersList[i]:
                        cluster_dict[j] = i
                modify_y = []
                for j in SCI.compRowIndexList + SCI.misRowIndexList:
                    modify_y.append(cluster_dict[j])
                origin_y, _ = SCI.modify_down_stream(cells)
                
                logging.debug("origin_y",origin_y)
                logging.debug("modify_y",modify_y)
                logging.debug("SCILN")
                if stop_check:
                    input("wait")
                
                self.purityAlg[algIndex] = Assist().purity(origin_y, modify_y)
                self.riAlg[algIndex] = Assist.RI(origin_y, modify_y)
                self.ariAlg[algIndex] = Assist.ARI(origin_y, modify_y)
                self.fmeasureAlg[algIndex],_,_ = Assist.f_measure(origin_y, modify_y)
                self.totalpurity[0][algIndex] += self.purityAlg[algIndex]
                self.totalri[0][algIndex] += self.riAlg[algIndex]
                self.totalari[0][algIndex] += self.ariAlg[algIndex]
                self.totalfmeasure[0][algIndex] += self.fmeasureAlg[algIndex]
            algTime = time() - startTime
            self.totalTime[0][algIndex] += algTime
        print("SCI_LN over")

    def SCI_LP(self):
        
        print("SCI_LP begin!")
        algIndex = 1
        if self.alg_flags[algIndex]:
            startTime = time()
            for repeat in range(self.RBEGIN, self.RBEGIN + self.REPEATLEN):
                seed = 100 + self.size * 1011 + repeat * 811
                cells, mask = self.dh.getMisCelMul(self.db, self.label_idx, self.size)

                SCI = SCILP(self.db, self.label_idx, self.exp_cluster, cells, mask)
                SCI.setK(self.K_SCI)
                SCI.setK_Candidate(self.K_Candidate_SCI)
                SCI.setEpsilon(self.epsilon_SCI)
                SCI.setDelCompRowIndexList([])
                SCI.mainSCI()
                cluster_dict = {}
                for i in range(len(SCI.clusterMembersList)):
                    for j in SCI.clusterMembersList[i]:
                        cluster_dict[j] = i
                modify_y = []
                for j in SCI.compRowIndexList + SCI.misRowIndexList:
                    modify_y.append(cluster_dict[j])
                origin_y, _ = SCI.modify_down_stream(cells)
                
                logging.debug("origin_y",origin_y)
                logging.debug("modify_y",modify_y)
                logging.debug("SCI-ILP")
                if stop_check:
                    input("wait")
                
                self.purityAlg[algIndex] = Assist().purity(origin_y, modify_y)
                self.riAlg[algIndex] = Assist.RI(origin_y, modify_y)
                self.ariAlg[algIndex] = Assist.ARI(origin_y, modify_y)
                self.fmeasureAlg[algIndex],_,_ = Assist.f_measure(origin_y, modify_y)
                self.totalpurity[0][algIndex] += self.purityAlg[algIndex]
                self.totalri[0][algIndex] += self.riAlg[algIndex]
                self.totalari[0][algIndex] += self.ariAlg[algIndex]
                self.totalfmeasure[0][algIndex] += self.fmeasureAlg[algIndex]
            algTime = time() - startTime
            self.totalTime[0][algIndex] += algTime
        print("SCI_LP over")

    def SCI_ILP(self):
        
        print(" SCI_ILP begin! ")
        algIndex = 2
        if self.alg_flags[algIndex]:
            startTime = time()
            for repeat in range(self.RBEGIN, self.RBEGIN + self.REPEATLEN):
                seed = 100 + self.size * 1011 + repeat * 811
                cells, mask = self.dh.getMisCelMul(self.db, self.label_idx, self.size)

                SCI = SCIILP(self.db, self.label_idx, self.exp_cluster, cells, mask)

                SCI.setK(self.K_SCI)
                SCI.setK_Candidate(self.K_Candidate_SCI)
                SCI.setEpsilon(self.epsilon_SCI)
                SCI.setDelCompRowIndexList([])
                SCI.mainSCI()

                self.totalTime[0][algIndex] += SCI.getAlgtime()

                cluster_dict = {}
                for i in range(len(SCI.clusterMembersList)):
                    for j in SCI.clusterMembersList[i]:
                        cluster_dict[j] = i
                modify_y = []
                for j in SCI.compRowIndexList + SCI.misRowIndexList:
                    modify_y.append(cluster_dict[j])
                origin_y, _ = SCI.modify_down_stream(cells)
                
                logging.debug("origin_y",origin_y)
                logging.debug("modify_y",modify_y)
                logging.debug("EXACT-SCI")
                if stop_check:
                    input("wait")
                
                self.purityAlg[algIndex] = Assist().purity(origin_y, modify_y)
                self.riAlg[algIndex] = Assist.RI(origin_y, modify_y)
                self.ariAlg[algIndex] = Assist.ARI(origin_y, modify_y)
                self.fmeasureAlg[algIndex],_,_ = Assist.f_measure(origin_y, modify_y)
                self.totalpurity[0][algIndex] += self.purityAlg[algIndex]
                self.totalri[0][algIndex] += self.riAlg[algIndex]
                self.totalari[0][algIndex] += self.ariAlg[algIndex]
                self.totalfmeasure[0][algIndex] += self.fmeasureAlg[algIndex]
            algTime = time() - startTime
            self.totalTime[0][algIndex] += algTime

        print("SCI_ILP over")

    def alg_exp(self):
        self.SCI_LN()
        self.SCI_LP()
        self.SCI_ILP()
        name1 = self.name + '_test'
        name2 = self.name
        columns = ["SCI_LN", "SCI_LP", "SCI_ILP"]
        purity_df = pd.DataFrame(( self.totalpurity / self.REPEATLEN), columns=columns)
        time_df = pd.DataFrame(( self.totalTime / self.REPEATLEN), columns=columns)
        ri_df = pd.DataFrame(( self.totalri / self.REPEATLEN), columns=columns)
        ari_df = pd.DataFrame(( self.totalari / self.REPEATLEN), columns=columns)
        fmeasure_df = pd.DataFrame(( self.totalfmeasure / self.REPEATLEN), columns=columns)
        if not os.path.exists(os.path.join("result/compare", name1)):
            os.makedirs(os.path.join("result/compare", name1))
            
            
        time_df.to_csv("result/" + "compare/" + name1 + "/" + name2 + "_time" + ".tsv", sep='\t',
                       float_format="%.5f", index=False)
           
        purity_df.to_csv("result/" + "compare/" + name1 + "/" + name2 + "_purity" + ".tsv", sep='\t',
                         float_format="%.5f", index=False)
            
        ri_df.to_csv("result/" + "compare/" + name1 + "/" + name2 + "_ri" + ".tsv", sep='\t',
                      float_format="%.5f",
                      index=False)
            
        ari_df.to_csv("result/" + "compare/" + name1 + "/" + name2 + "_ari" + ".tsv", sep='\t',
                      float_format="%.5f", index=False)
            
        fmeasure_df.to_csv("result/" + "compare/" + name1 + "/" + name2 + "_f1" + ".tsv", sep='\t',
                           float_format="%.5f", index=False)


        print("all over!")


if __name__ == '__main__':
    ct = CompTest("../data/crx/crx.data", label_idx=15, exp_cluster=2, header=None, index_col=None)
    ct.alg_exp()
    
