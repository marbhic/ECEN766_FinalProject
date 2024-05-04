# This file filters the test dataset to only breast cancer cell lines. Breast Cancer cell lines are given in the breast_cancer_lines.csv

import pandas as pd

breast_lines = pd.read_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1d\\breast_cancer_lines.csv")   # Breast cancer Cell Lines
test_set = pd.read_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1c\\Fig1c_Precily_pathways\\test_data_attempt2.csv") # Location of Test dataset that should be filtered

breast_lines_list = breast_lines["stripped_cell_line_name"].to_list()
breast_cancer_df = test_set[test_set['stripped_cell_line_name'].isin(breast_lines_list)]
breast_cancer_df.to_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1d\\breast_cancer_testset.csv", index = False) 
