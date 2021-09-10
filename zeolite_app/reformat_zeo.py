#!/usr/bin/env python
from load_data import GetZeoliteTsv
import pandas as pd
import pickle

class FormatZeo(GetZeoliteTsv):

    def __init__(self, new_entries, template = "zeolite_template.tsv"):

        super().__init__(template)
        self.new_entries = new_entries
        self.rows, self.cols = self.zeolite_df.shape
        self.set_dtypes()
        self.zeolite_df = pd.concat([self.zeolite_df, self.new_entries])
        self.encode_categorical()       

    def get_entry(self):
        
        return self.zeolite_df.iloc[self.rows:,:]

    def save_entry(self, outfile):
        
        self.zeolite_df = self.zeolite_df.iloc[self.rows:,:]
        self.save_zeo(outfile)
    
            
# zeo  = FormatZeo( pd.read_csv("1.csv", sep = "\t"))
# print(zeo.get_entry())
# zeo.save_entry("Entry.tsv")
