#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=Panc-AdenoCA
#python3 3_Intersect_paz_cancertypes.py --cancer_types=Panc-AdenoCA --decrement_by=23
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=Panc-AdenoCA

#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=Panc-AdenoCA
#python3 3_Intersect_paz_cancertypes.py --cancer_types=Panc-AdenoCA --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=Panc-AdenoCA

#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=Lung-SCC
#python3 3_Intersect_paz_cancertypes.py --cancer_types=Lung-SCC --increment_by=5 --max_samples=45
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=Lung-SCC
#Rscript align_mutations_to_ranges.R --cancer_type=Lung-SCC


#cancer_type="Lung-AdenoCA"
#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types $(cancer_type)
#python3 3_Intersect_paz_cancertypes.py --cancer_types $cancer_type --decrement_by=3
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types $cancer_type

#cancer_type="Liver-HCC"
#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types $cancer_type
#python3 3_Intersect_paz_cancertypes.py --cancer_types $cancer_type --increment_by=5 --max_samples=45
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types $cancer_type
#Rscript align_mutations_to_ranges.R --cancer_type=$cancer_type

cancer_type="CNS-GBM"
python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types $cancer_type
python3 3_Intersect_paz_cancertypes.py --cancer_types $cancer_type --increment_by=5 --max_samples=35
python3 4_AssembleCout_paz_Cancergroup.py --cancer_types $cancer_type
Rscript align_mutations_to_ranges.R --cancer_type=$cancer_type


#cancer_type="Kidney-ChRCC"
#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types $cancer_type
#python3 3_Intersect_paz_cancertypes.py --cancer_types $cancer_type --increment_by=5 --max_samples=40
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types $cancer_type
#Rscript align_mutations_to_ranges.R --cancer_type=$cancer_type


#cancer_type="ccRCC"
#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types $cancer_type
#python3 3_Intersect_paz_cancertypes.py --cancer_types $cancer_type --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types $cancer_type

#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=Skin-Melanoma
#python3 3_Intersect_paz_cancertypes.py --cancer_types=Skin-Melanoma --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=Skin-Melanoma
#Rscript align_mutations_to_ranges.R --cancer_type=Skin-Melanoma


#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=ColoRect-AdenoCA
#python3 3_Intersect_paz_cancertypes.py --cancer_types=ColoRect-AdenoCA --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=ColoRect-AdenoCA

#cancer_type="CNS-Medullo"
#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=CNS-Medullo
#python3 3_Intersect_paz_cancertypes.py --cancer_types=CNS-Medullo --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=CNS-Medullo

#python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types=CNS-PiloAstro
#python3 3_Intersect_paz_cancertypes.py --cancer_types=CNS-PiloAstro --increment_by=5 --max_samples=50
#python3 4_AssembleCout_paz_Cancergroup.py --cancer_types=CNS-PiloAstro
