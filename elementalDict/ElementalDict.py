#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


database_electronic_properties = pd.read_excel("database_electronic properties.xlsx",  header = 1, index_col = 0)


# In[3]:


elemental_property = pd.read_excel("oliynyk-elemental-property-list.xlsx",  header = 0, index_col = 0)


# In[4]:


elemental_property.columns


# In[5]:


def get_weighted_atomic_radius(row):

    print(row['Alloy'])
    
    clean_row = row.dropna()
    
    #get only element in the reference
    valid_indexes = elemental_property.index.intersection(clean_row.index)

    #Temporay data.frame
    df = pd.DataFrame()
    df['atomic_fraction'] = row.loc[valid_indexes]
    df['atomic_radius'] = elemental_property.loc[valid_indexes, ['Atomic\nradius calculated']]
    df['weighted_atomic_radius'] = df.atomic_fraction * df.atomic_radius
    print(df) 
    return df['weighted_atomic_radius'].sum()


# In[6]:


database_electronic_properties['weighted_atomic_radius'] = database_electronic_properties.apply(get_weighted_atomic_radius, axis=1)


# In[7]:


database_electronic_properties[['Alloy','Precusor','weighted_atomic_radius']]


# In[ ]:




