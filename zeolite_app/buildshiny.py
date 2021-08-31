#!/usr/bin/env python
from load_data import GetZeoliteTsv
import pandas as pd
import argparse



class Table2Buttons(GetZeoliteTsv):

    def __init__(self, zeolite_filename):
        super().__init__(zeolite_filename)
        self.set_dtypes()
        self.metals = ['Na+','Ag+', 'Cu+', 'Ce+4', 'Cs+2', 'Ni+2']
        self.metal_vars = [ col  for col in self.zeolite_df.columns if col.startswith("x_") or col.startswith("R_")]
        self.metal_props = dict([ (var, self.zeolite_df[var].max()) for var in self.metal_vars ])
        self.selectInput = [ col for col,dtype  in dict(self.zeolite_df.dtypes).items() if dtype == 'O' ]
        self.sliderInput =  [ col for col,dtype  in dict(self.zeolite_df.dtypes).items() if dtype == 'float64' ]
        
        
       
    def get_Input(self):

        for input_col1 in self.selectInput:
            choices = self.zeolite_df[input_col1].unique()
            self.tag_selectInput(input_col1, tuple(choices))

        self.tag_selectInput("Metals", tuple(self.metals), multiple ='TRUE')    
            
        for input_col2 in self.sliderInput:
            metrics = self.zeolite_df[input_col2].describe()
            if input_col2 in  self.metal_vars:continue
            if input_col2 == 'Capacity':continue
            
            self.tag_sliderInput(input_col2, metrics['min'], metrics['max'], metrics['mean'])

    
    def tag_selectInput(self, inputId, choices, multiple ='FALSE'):
        
        print("""selectInput("{inputId}", 
"{inputId}:", 
choices = c{choices},
multiple = {multiple}),\n""".format(**{"inputId": inputId.title(), "choices": repr(choices), "multiple":multiple}))

        


    def tag_sliderInput(self, inputId, vmin, vmax, mean):

        print("""sliderInput("{inputId}",
"{inputId}",
{min},
{max},
{mean}),\n""".format(**{"inputId": inputId.title(),"min": vmin,"max": vmax, "mean": mean})) 




        
        

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="build shinu app buttons")
    parser.add_argument('zeolite_file',help = "zeolite file")
    args = parser.parse_args()
    t2b = Table2Buttons(args.zeolite_file)
    t2b.set_dtypes()
    t2b.get_Input()
    #for metal in ['Na+','Ag+', 'Cu+', 'Ce+4', 'Cs+2', 'Ni+2', 'x_Na+', 'x_Ag+', 'x_Cu+','x_Ce+4', 'x_Cs+2', 'x_Ni+2', 'R_Na+', 'R_Ag+', 'R_Cu+', 'R_Ce+4','R_Cs+2', 'R_Ni+2']:
    #     print(meta)
   
    #t2b.get_selectInput()
    #t2b.get_sliderInput()
    
