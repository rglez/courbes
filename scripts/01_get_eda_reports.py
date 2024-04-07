# Created by roy.gonzalez-aleman at 07/04/2024
"""
Script dealing with the exploratory data analysis report
"""

import pandas as pd
from ydata_profiling import ProfileReport, compare

descriptor_file_a = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_WT/inter/Roll.txt'
descriptor_file_b = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_SNO/inter/Roll.txt'
descriptor_file_c = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_SOH/inter/Roll.txt'
df_a = pd.read_table(descriptor_file_a, header=0, sep='\s+')
df_b = pd.read_table(descriptor_file_b, header=0, sep='\s+')
df_c = pd.read_table(descriptor_file_c, header=0, sep='\s+')

# prof_a = ProfileReport(df_a, title='WT')
prof_b = ProfileReport(df_b, title='SNO')
prof_c = ProfileReport(df_b, title='SOH')
# prof.to_file(output_file='output.html')

comparison_report = compare([prof_b, prof_c])
comparison_report.to_file("comparison.html")
