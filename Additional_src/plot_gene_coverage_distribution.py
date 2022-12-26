#!/usr/bin/env python3
# coding=utf-8

"""
Этот скрипт нужен, чтобы делать синаплоты покрытия генов однокопийных и многокопийных ортогрупп.

На вход он берёт 
1) Файл logs.txt , который делает calculate_AG.py
2) Размер точки синаплота. Нормальный размер 2.5, но если ортогрупп очень много, то можно и уменьшить.
3) Какое взять максимальное покрытие. Или слово "auto" или число. Если пользователь указал auto, то скрипт сам выбирает максимальное покрытие для диаграммы (то есть, максимальное значение, которое будет отображено на оси Y), а если пользователь указал число, то скрипт использует такое покрытие.
4) Путь к выходному файлу с картинкой без расширения. 

Скрипт создаёт два файла с одной и той же картинкой — один в формате svg и один в формате png.

Пример использования:
python3 plot_gene_coverage_distribution.py ./Soldierfish_assembly_results/AG_calculation_for_max_divergence_0.216452/logs.txt 2.5 auto ./Soldierfish_assembly_results/AG_calculation_for_max_divergence_0.216452/gene_coverage_distribution

Соответственно, в папке ./Soldierfish_assembly_results/AG_calculation_for_max_divergence_0.216452/ создадутся файлы gene_coverage_distribution.svg и gene_coverage_distribution.png .
"""

import sys
import os
import subprocess
import re
import statistics
import pandas
import plotnine

f_infile = open(sys.argv[1], "r")
n_point_size = float(sys.argv[2])
s_maximum_coverage_to_draw = sys.argv[3]
s_path_to_the_output_plot_without_extension = sys.argv[4]
s_path_to_the_output_svg_file = s_path_to_the_output_plot_without_extension + ".svg"
s_path_to_the_output_png_file = s_path_to_the_output_plot_without_extension + ".png"

#составлю список с покрытиями генов и список такой же длины, в котором фразы "Single-copy orthogroups" и "Multicopy orthogroups" — по одной фразе для каждого значения покрытия. Эти два списка имеют одинаковую длину.

l_orthogroup_copy_numbers = [] #Значения могут быть "Single-copy orthogroups" или "Multicopy orthogroups"
l_coverages = [] #медианные покрытия экзонов генов.

for s_line in f_infile:
	#For the gene contig_1623:64878-86622 from a single-copy orthogroup 72550at7898, the median coverage is 54.0. It was calculated using 1002 positions.
	o_regular_expression_results = re.search(r"^For the gene .+ from a (.+) orthogroup .+ the median coverage is ([\d\.\-\+eE]+)\.\s", s_line)
	
	if o_regular_expression_results:
		s_orthogroup_copy_number = o_regular_expression_results.group(1)
		n_coverage = float(o_regular_expression_results.group(2))
		
		if (s_orthogroup_copy_number == "single-copy") or (s_orthogroup_copy_number == "single copy"): #в Mabs 2.9 было "single copy", в Mabs 2.10 стало "single-copy".
			l_orthogroup_copy_numbers.append("Genes from\nsingle-copy\northogroups")
		if s_orthogroup_copy_number == "multicopy":
			l_orthogroup_copy_numbers.append("Genes from\nmulticopy\northogroups")
		
		l_coverages.append(n_coverage)

f_infile.close()

#print(l_orthogroup_copy_numbers)
#print(l_coverages)

#конвертирую список в data frame
o_data_frame = pandas.DataFrame(data = {"": l_orthogroup_copy_numbers, "Sequencing coverage": l_coverages})

#print(o_data_frame)

#Если пользователь указал s_maximum_coverage_to_draw = "auto", то диапазон по оси ординат ставлю 0, max(2 * median(y), 1.05 * max(y)). Такое максимальное значение я ставлю потому, что мне хочется (для красоты), чтобы центр распределения был по середине графика. Но, вместе с тем, если самая высокая точка выше, чем две медианы (2 * median(y)), то я сделаю максимум ещё выше: на 5% выше, чем самая высокая точка. 5% добавляю, чтобы график несколько лучше выглядел.
if s_maximum_coverage_to_draw == "auto":
	n_maximum_coverage_to_draw = max(2 * statistics.median(l_coverages), 1.05 * max(l_coverages))
#если пользователь указал максимальное покрытие числом.
else:
	n_maximum_coverage_to_draw = float(s_maximum_coverage_to_draw)

#цвет — экселевский синий. Д
#plotnine.scales.scale_y_continuous(expand = (0, 0)) я добавляю, чтобы не было отступа снизу по оси ординат.
#text = plotnine.element_text(linespacing = 1.1) я добавляю потому, что иначе расстояние между строками текста слишком маленькое.

o_plot = plotnine.ggplot(o_data_frame) + plotnine.aes(x = "", y = "Sequencing coverage") + plotnine.geom_sina(color = "#4472C4", scale = "area", size = n_point_size) + plotnine.scales.scale_x_discrete(limits = ["Genes from\nsingle-copy\northogroups", "Genes from\nmulticopy\northogroups"]) + plotnine.scales.scale_y_continuous(expand = (0, 0), limits = (0, n_maximum_coverage_to_draw)) + plotnine.theme_light() + plotnine.theme(text = plotnine.element_text(linespacing = 1.1))

#Делаю картинки
o_plot.save(filename = s_path_to_the_output_svg_file, format = "svg", height = 10, width = 10, units = 'cm')
o_plot.save(filename = s_path_to_the_output_png_file, format = "png", height = 10, width = 10, units = 'cm', dpi = 300)



