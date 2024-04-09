import pandas as pd
import numpy as np


# original dataset
old_dataset = pd.read_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1d\\Test_Set\\Test_set.csv")
# cell line names for mutations
new_data = pd.read_csv("C:\\Users\\marbhic\\Downloads\\sample_info.csv")
# mutation data
mutation_data = pd.read_csv("C:\\Users\\marbhic\\Downloads\\filtered_mutations.csv")

# merge cell line names and mutation data
df_merged = mutation_data.merge(new_data[['DepMap_ID', 'stripped_cell_line_name']])
cell_lines = df_merged['stripped_cell_line_name'].tolist()

# merge mutation data and original dataset
x = old_dataset.copy()
x.rename(columns = {'0':'stripped_cell_line_name'}, inplace = True) 

# type casting data to decrease size of file 
x[x.select_dtypes(np.float64).columns] = x.select_dtypes(np.float64).astype(np.float32)
df_merged[df_merged.select_dtypes(np.float64).columns] = df_merged.select_dtypes(np.float64).astype(np.float32)
df_merged[df_merged.select_dtypes(np.int64).columns] = df_merged.select_dtypes(np.int64).astype(np.int32)
merge_mutation = x.merge(df_merged, on= 'stripped_cell_line_name')
merge_mutation.to_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1c\\Fig1c_Precily_pathways\\test_data_attempt2.csv")


