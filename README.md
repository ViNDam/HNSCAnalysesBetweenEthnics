# HNSCMethylationBetweenEthnics
This code is a part of the published paper related to immune landscapes in Afican Americanas versus Whites.

Chaudhary, S., Dam, V., Ganguly, K., Sharma, S., Atri, P., Chirravuri-Venkata, R., ... & Batra, S. K. (2020). Differential mutation spectrum and immune landscape in African Americans versus Whites: A possible determinant to health disparity in head and neck cancer. Cancer letters, 492, 44-53.

This code download TCGA Methylation data from manifest file, generate analyses including Beta Density Plot, Mean Values By Condition Plot, and adjusted Wilcox p value between "Negative_black_or_african_american", "Negative_white" and "Positive_white". 

Then, using adjusted wilcox p_value and fold change difference between 2 ethnic groups, it classifies the genes into three groups "Fold Change" (FC > 1), :SignificantFC" (p<0.05 & FC > 1) and "Significant (p<0.05). 

This code helps identify the differences in methylation and expression between ethnic groups in head and neck cancer. The user can apply similar analyses for other cancer types. 

# CODE

Methylation.R

The output table includes cg ID, gene symbol, p value, log2 fold change value , and significant status will be created for each group. 
Output plots such as significant genes plot , survival analysis plot, enrichment pathway plot, Beta Density Plot, Mean Values By Condition Plot will be produced. 

# RESULT of HSNC cancer
<img src="outputs/HSNC.png" width="400" height="400">

<img src="outputs/BetaDensity.png" width="400" height="400">

<img src="outputs/survival_pos_vs_neg_W.png" width="400" height="400">

<img src="outputs/Enrichment pathway PW.png" width="400" height="400">

<img src="outputs/Mutational_burden.png" width="400" height="400">



