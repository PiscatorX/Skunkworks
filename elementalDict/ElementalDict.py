#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


pd.set_option('display.max_rows', None)


# In[3]:


database_electronic_properties = pd.read_excel("database_electronic properties.xlsx",  header = 1, index_col = 0)


# In[4]:


elemental_property = pd.read_excel("oliynyk-elemental-property-list.xlsx",  header = 0, index_col = 0)


# In[5]:


elemental_property.columns


# In[6]:


def get_weighted_atomic_radius(row):

    alloy = (row)
    clean_row = row.dropna()
    
    #get only element in the reference
    valid_indexes = elemental_property.index.intersection(clean_row.index)
    
    atomic_radius = elemental_property.loc[valid_indexes, ['Atomic\nradius calculated']]
    
    return atomic_radius['Atomic\nradius calculated'].sum()


# In[7]:


database_electronic_properties['weighted_atomic_radius'] = database_electronic_properties.apply(get_weighted_atomic_radius, axis=1)


# In[8]:


database_electronic_properties[['Alloy','Precusor','weighted_atomic_radius']]


# In[ ]:




