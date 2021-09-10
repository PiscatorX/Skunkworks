#!/usr/bin/env python
from sklearn import metrics
import pandas as pd
import numpy as np
import pickle
import joblib


class ZeoRandomForest(object):

    def __init__(self,
                 entries_df,
                 StandardScaler = "project-rodger-dogder-StandardScaler.obj",
                 RF_model  = "project-rodger-dogder-RF_model.obj",
                 Xtest = "X_test.pickle",
                 ytest = "y_test.pickle"):

        self.RF_regressor_model = joblib.load(RF_model)
        self.RF_params = self.RF_regressor_model.get_params()
        self.sc = joblib.load(StandardScaler)
        self.Xtest = pd.read_pickle(Xtest)
        self.Xtest.to_csv("SampleX.tsv", sep = "\t")
        entries_df.columns = self.Xtest.columns
        self.X_real = self.sc.transform(entries_df)
        self.Xtest = self.sc.transform(pd.read_pickle(Xtest))
        self.ytest = pd.read_pickle(ytest)
        
    
    def RF_predict(self):

        return self.RF_regressor_model.predict(self.X_real)

    
    
    def check_r2(self):
        
         y_pred = self.RF_regressor_model.predict(self.Xtest)
         
         return metrics.r2_score(self.ytest, y_pred)


# if  __name__ ==  '__main__':
#      entries = pd.read_pickle("X_test.pickle")
#      c1 = set(entries.columns)
#      entries1 = pd.read_csv('1.csv')
#      c2 = set(entries.columns)
     
#      zeo = ZeoRandomForest(entries)
#      print(zeo.RF_predict())
#      #print(zeo.check_r2())
