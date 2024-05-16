# Supplementary code repository for: Unravelling the transcriptome of the human tuberculosis lesion and its clinical implications.

# Authors

Kaori L. Fonseca , Juan José Lozano, Albert Despuig, Dominic Habgood-Coote, Julia Sidorova, Lilibeth Arias, Álvaro Del Río-Álvarez, Juan Carrillo-Reixach,Aaron Goff, Leticia Muraro Wildner, Shota Gogishvili, Keti Nikolaishvili, Natalia Shubladze, Zaza Avaliani, Pere-Joan Cardona, Federico Martinón-Torres,Antonio Salas,Alberto Gómez-Carballa, Carolina Armengol, Simon J Waddell, Myrsini Kaforou, Anne O'Garra, Sergo Vashakidze, Cristina Vilaplana

#  Introduction:

Three separate R scripts and corresponding input files are dedicated to the main analysis steps to reproduce the most important published results. 

#  Analysis workflow R - scripts

* script1.R - Required for figures 2a, 2b and 3b. DEG determination and compartment enrichment analysis.
* script2.R - Required for WGCNA processed generated modules.
* script3.R - Required for Qusage enrichment from generated modules 


# Repository files:

* TBrnaseq.rds  - RDS includind feature counts data and minimal annotation used por 
* phenotype.xls - metadata including the cellular compartment and some clinical variables. 
* GranulomeModules.gmt - Derived 21 gene modules related to TB compartments (derived from script2.R) 

 
