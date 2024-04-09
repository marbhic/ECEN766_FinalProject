# Predicting Drug Response in Breast Cancer Using Deep Learning Neural Networks

<h2>Overview</h2>
<h2>Code Overview</h2>
<p>This GitHub repository will provide all the required data and code to run the modified version of Precily. The following is the breakdown of how the project comes together</p>
<ol>
  <li><b>1.MutationDataMerge</b></li>
  <p>This folder provides all the data and code to create the new test and trained dataset. Within the folder are the original dataset, the mutation data, the code to create the new training and test set. Mutation data merge with pathway scores and drug descriptors was done in 3 steps. <br> </p>
  <b>1.MutationDataMerge/1.Mutation_DataFiltering.R</b> <br>
  This file takes mutation data provided by the DepMap portal in their 20Q1 release (<b>CCLE_mutations.csv</b>), trims unnecessary columns, and overlaps mutations with Cosmic Cancer Gene Cencus data (<b>cancer_gene_census.csv</b>) where it was first converted to GRCh38 format using the liftOver Library (<b>hg19ToHg38.over.chain</b>) and generates a table of size 1580 rows by 495 columns where each row indicated the cell line and indications if a certain mutation was present or not (<b>FILTERD_DepMap_21Q2_Mutations_by_Cell.csv</b>). <br><br>
  <b>1.MutationDataMerge/2.DataSetMerge.py</b><br>
  This file merges the mutation matrix created in 1.Mutation_DataFiltering.R with the original Pathway+DrugDescriptors   
  training and test set (<b>/Orignal_dataset</b>). <br><br>
  <b>Manual Update</b><br>
  The merged training and test set created in part 2 required additional manual configuration, the IC50 column was moved to   
  the last coloumn. The completely merged training and test set are located at <b>1.MutationDataMerge/Mutation+OldData/</b>
  <li>Precily Training</li>
  <li>Evaluation</li>
</ol>
