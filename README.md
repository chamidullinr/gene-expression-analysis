# Gene Expression Analysis
Assignment on learning from Gene Expression Data and applying Statistical Analysis methods.

## Introduction
Statistical Analysis is an exiting field that can be considered as a cousin of Machine Learning.
While Machine Learning is often placed under Computer Science, Statistical Analysis is closer to Mathematics.
Both approaches focus on finding interpretation and relationships in the data 
(tasks such as data mining and modeling), but these fields differ in the workflow.
There is not a clear line that separates domains of Statistical Analysis and Machine Learning 
and there is a lot of algorithms that overlap in both approaches.

Usual practice in Machine Learning is to pre-process the data (i.e. apply cleaning, features scaling, etc.),
then employ Machine Learning algorithms (for tasks such as classification, regression, clustering, etc.),
and in the end evaluate performance of applied algorithms on the given data.
In other words, in Machine Learning one first applies algorithms and then reasons about the meaning of the results - 
whether the applied algorithm was appropriate for given problem and what can be improved.

In contrast, Statistical Analysis tackles the problem from the different perspective.
First, one makes assumptions about the data and then applies the relevant algorithms.
In this approach, there are often utilized methods such as hypothesis testing and correlation analysis.

To myself personally, Statistical Analysis is an interesting field where one needs to have a good understanding 
of methods, statistical tests and predictive models.
Often such methods cannot be just applied as a black-box, it requires testing 
whether it is possible to employ the method in the first place.
And also it requires understanding the data and its distribution.

This is not to undermine Machine Learning. Both Machine Learning and Statistical Analysis are interesting fields
with different focuses and challenges. Just sometimes, in Machine Learning, some aspects of data analysis and modeling are neglected. 

In this assignment, the goal is to create hypotheses about the data and relationships between variables. 
And then verify them with statistical methods. The task focuses on the data from Bioinformatics.

## Task
The provided dataset covers gene expression data that contain expression measurements 
of a large pool of genes over a reasonable sample of individuals. 
The individuals are split into several well-defined groups such as:
* healthy controls,
* stationary patients,
* chronically diseased patients and
* recovered patients.
  
The goal is to learn from the gene expression data and understand relationships 
among activity of the individual genes, and the health condition of an individual.

## Solution
The completed task is located in the Jupyter Notebook **gene_expression.ipynb**.
The report contains several sections that cover the following tasks:
1) Exploratory analysis of the gene expression dataset
    * correlation analysis
    * feature scaling
    * dimensionality reduction
    * clustering
2) Differential expression
    * ANOVA and Tukey HSD as post-hoc ANOVA
3) Phenotype predictive model
    * Linear Discriminant Analysis (LDA)
    * Logistic Regression
    * Random Forest
    * Support Vector Machine (SVM) - with a linear kernel

## Package Requirements
General packages for working with data:
* numpy
* pandas

Scientific packages for Machine Learning and Statistical Analysis:
* scipy
* scikit-learn
* statsmodels

Visualization packages:
* matplotlib
* seaborn


## Authors
**Rail Chamidullin** - chamidullinr@gmail.com  - [Github account](https://github.com/chamidullinr)
