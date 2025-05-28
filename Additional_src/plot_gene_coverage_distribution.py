#!/usr/bin/env python3
# coding=utf-8

"""
plot_gene_coverage_distribution.py is a part of the genome assembly suite Mabs. See https://github.com/shelkmike/Mabs

Notes:
1) All comments, except this one, are in Russian. Sorry, but it's somewhat easier for me to write in Russian than in English. To understand some comment, you can use Google Translate. Names of variables are usually self-explanatory, so it is often possible to understand the meaning of a piece of code without comments. In case of trouble understanding code, ask a question at https://github.com/shelkmike/Mabs/issues .
2) Throughout the code I use a Hungarian notation, which means I denote meaning of a word by using special prefixes. In particular:
s_ - variables that contain strings
n_ - variables that contain numbers
l_ - lists
d_ - dictionaries
f_ - file handlers
o_ - more complex objects
Nested data structures are denoted by several letters. For example, dl_ are dictionaries of lists and ll_ are lists of lists.

For details on how to use this script, run "python3 plot_gene_coverage_distribution.py --help".

Этот скрипт нужен, чтобы делать синаплоты покрытия генов однокопийных и многокопийных ортогрупп. Каждая точка это ген. Точки синим цветом.
Кроме этого, он может делать синаплоты для самих ортогрупп. Каждая точка это ортогруппа, а значение — среднее покрытие её генов. Точки зелёным цветом.
На вход скрипт берёт папку, выдаваемую calculate_AG. Подробнее можно посмотреть, запустив скрипт с ключом --help.
"""

import sys
import os
import subprocess
import re
import statistics
import pandas
import plotnine

#делаю парсинг аргументов командной строки. Можно было бы использовать argparse, но когда я делаю это без библиотек, то больше возможностей для того, чтобы сделать интерфейс таким, какой мне нравится.

s_command_line = " ".join(sys.argv) #команда, которой запущен скрипт, в одну строку.
s_command_line_reduced = s_command_line #то же, что s_command_line, но после того, как я распаршу какой-нибудь аргумент, я удалю его из этой строки. Если останется какой-то нераспарсенный аргумент, значит пользователь ввёл неизвестные скрипту аргументы, и нужно выдать ошибку.

#инициализирую исходные значения переменных
s_path_to_the_folder_with_calculate_AG_results = ""
s_path_to_the_output_folder = "Coverage_distribution"
s_maximum_coverage_to_draw = "auto"
n_point_size = 2.5
s_distribution_area_type = "equal"
s_make_orthogroup_coverage_plots = "false" 
n_number_of_bins_in_the_sina_plots = 50
n_random_seed = 1

l_errors_in_command_line = [] #список ошибок в командной строке. Если пользователь совершил много ошибок, то я напишу на экран сразу все, а не только первую встреченную.

#если нет ни одного аргумента командной строки, или есть аргумент командной строки --help, то печатаю хелп
if (len(sys.argv) == 1) or re.search(r"\s\-\-help", s_command_line):
	print("""This script plots per-gene and per-orthogroup coverage distributions as sina plots (https://en.wikipedia.org/wiki/Sina_plot). \nAs input, it takes a folder produced by calculate_AG. For details, see https://github.com/shelkmike/Mabs.

Main options:
1) --calculate_AG_folder       Path to a folder created by calculate_AG. This is the only mandatory option.
2) --output_folder        Path to the output folder. The default value is "Coverage_distribution".
3) --maximum_coverage_to_draw     Maximum value for the vertical axis. Should be either a number or the word "auto". In the latter case the maximum value is determined automatically based on gene coverage distribution. The default value is "auto".
4) --point_size     Point size for the sina plots, in millimeters. The default value is 2.5.
5) --distribution_area_type     Whether area of distributions of single-copy and multicopy orthogroups should be equal to each other (the value is "equal"), or proportional to the number of points (the value is "proportional"). The default value is "equal".
6) --number_of_bins     Number of bins in the sina plots. Should be an integer. The default value is 50.
7) --seed     A pseudorandom number generator's seed that affects jitter of points in the sina plots. Should be an integer. Usually this parameter shouldn't be changed. The default value is 1.
8) --make_orthogroup_coverage_plots     When --make_orthogroup_coverage_plots is set to "true", in addition to gene sina plots this script will also make orthogroup sina plots. In gene sina plots every point is a gene, while in orthogroup sina plots every point is an orthogroup. For example, if a BUSCO orthogroup has 2 genes, in the orthogroup sina plot they will be represented by 1 point that has coverage equal to the arithmetic mean of coverages of these two genes. The orthogroup sina plots will be in files orthogroup_coverage_distribution.png and orthogroup_coverage_distribution.svg. The value of this option should be either "true" or "false". The default value of this option is "false".

Informational options:
9) --help        Print this help.

Example:
plot_gene_coverage_distribution.py --calculate_AG_folder AG_calculation_results
	""")
	sys.exit()

#смотрю, указал ли пользователь в командной строке путь к папке с результатами calculate_AG.
o_regular_expression_results = re.search(r" --calculate_AG_folder (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_path_to_the_folder_with_calculate_AG_results = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь в командной строке путь к выходной папке.
o_regular_expression_results = re.search(r" --output_folder (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_path_to_the_output_folder = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь в командной строке максимальное покрытие, которое нужно нарисовать.
o_regular_expression_results = re.search(r" --maximum_coverage_to_draw ([\d\.eE\-\+]+|auto)", s_command_line_reduced)
if o_regular_expression_results:
	s_maximum_coverage_to_draw = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь в командной строке размер точки для рисунков.
o_regular_expression_results = re.search(r" --point_size ([\d\.eE\-\+]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_point_size = float(o_regular_expression_results.group(1))
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь в командной строке тип распределения (с равной площадью или с разной площадью).
o_regular_expression_results = re.search(r" --distribution_area_type (equal|proportional)", s_command_line_reduced)
if o_regular_expression_results:
	s_distribution_area_type = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь в командной строке количество бинов в синаплоте.
o_regular_expression_results = re.search(r" --number_of_bins ([\d\.eE\+]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_number_of_bins_in_the_sina_plots = int(o_regular_expression_results.group(1))
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь сид генератора случайных чисел, регулирующего положение точек в синаплоте.
o_regular_expression_results = re.search(r" --seed ([\d\.eE\+]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_random_seed = int(o_regular_expression_results.group(1))
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, использовал ли пользователь --make_orthogroup_coverage_plots.
o_regular_expression_results = re.search(r" --make_orthogroup_coverage_plots (true|false|True|False|TRUE|FALSE)", s_command_line_reduced)
if o_regular_expression_results:
	s_make_orthogroup_coverage_plots = o_regular_expression_results.group(1)
	#перевожу в нижний реестр.
	s_make_orthogroup_coverage_plots = s_make_orthogroup_coverage_plots.lower()
					
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#проверяю, не ввёл ли пользователь какие-то несуществующие опции. Это я определяю по тому, что после того, как я распарсил все команды, в строке s_command_line_reduced осталось что-то, кроме названия исполняемого файла скрипта.
s_command_line_reduced = re.sub(r"^.*?plot_gene_coverage_distribution(\.py)?\s*", "", s_command_line_reduced)
if s_command_line_reduced != "":
	l_errors_in_command_line.append("You have provided some options which the script doesn't know: " + s_command_line_reduced)

#Проверяет, есть ли папка с результатами calculate_AG и есть ли в ней файл log.txt. Если папки нет (или папка есть, но в ней нет файла), то это ошибка.
if not os.path.exists(s_path_to_the_folder_with_calculate_AG_results):
	l_errors_in_command_line.append("The input folder does not exist.")

if os.path.exists(s_path_to_the_folder_with_calculate_AG_results) and (not os.path.exists(s_path_to_the_folder_with_calculate_AG_results + "/log.txt")):
	l_errors_in_command_line.append("The file log.txt does not exist the input folder.")

#Пишу список ошибок в командной строке.
#Если ошибка была всего одна.
if len(l_errors_in_command_line) == 1:
	print("There was an error in the command line:")
	print(l_errors_in_command_line[0])
	sys.exit()
#Если было больше одной ошибки.
if len(l_errors_in_command_line) > 1:
	print("There were errors in the command line:")
	n_error_number = 0 #порядковый номер ошибки. Считается от 1.
	for s_error_text in l_errors_in_command_line:
		n_error_number += 1
		print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
	sys.exit()

#Печатаю пустую строку, как разделитель
print("")

#Проверяю, существует ли выходная папка. Если её не существует, то создаю её.
if not os.path.exists(s_path_to_the_output_folder):
	os.mkdir(s_path_to_the_output_folder)	


#Сначала сделаю графики для генов

#составлю список с покрытиями генов и список такой же длины, в котором фразы "Single-copy orthogroups" и "Multicopy orthogroups" — по одной фразе для каждого значения покрытия. Эти два списка имеют одинаковую длину.
f_infile = open(s_path_to_the_folder_with_calculate_AG_results + "/log.txt", "r")
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

#цвет — экселевский синий.
#plotnine.scales.scale_y_continuous(expand = (0, 0)) я добавляю, чтобы не было отступа снизу по оси ординат.
#text = plotnine.element_text(linespacing = 1.1) я добавляю потому, что иначе расстояние между строками текста слишком маленькое.

if s_distribution_area_type == "equal":
	o_plot = plotnine.ggplot(o_data_frame) + plotnine.aes(x = "", y = "Sequencing coverage") + plotnine.geom_sina(color = "#4472C4", scale = "area", size = n_point_size, random_state = n_random_seed, bins = n_number_of_bins_in_the_sina_plots) + plotnine.scales.scale_x_discrete(limits = ["Genes from\nsingle-copy\northogroups", "Genes from\nmulticopy\northogroups"]) + plotnine.scales.scale_y_continuous(expand = (0, 0), limits = (0, n_maximum_coverage_to_draw)) + plotnine.theme_light() + plotnine.theme(text = plotnine.element_text(linespacing = 1.1))
if s_distribution_area_type == "proportional":
	o_plot = plotnine.ggplot(o_data_frame) + plotnine.aes(x = "", y = "Sequencing coverage") + plotnine.geom_sina(color = "#4472C4", scale = "count", size = n_point_size, random_state = n_random_seed, bins = n_number_of_bins_in_the_sina_plots) + plotnine.scales.scale_x_discrete(limits = ["Genes from\nsingle-copy\northogroups", "Genes from\nmulticopy\northogroups"]) + plotnine.scales.scale_y_continuous(expand = (0, 0), limits = (0, n_maximum_coverage_to_draw)) + plotnine.theme_light() + plotnine.theme(text = plotnine.element_text(linespacing = 1.1))
	
#Делаю картинки.
#"verbose = False" использую, чтобы Plotnine не писал warnings про то, что он сохраняет картинки. Почему-то когда Plotnine сохраняет картинки он всегда об этом пишет, используя слово "Warning". Это может смутить пользователя.
o_plot.save(filename = s_path_to_the_output_folder + "/gene_coverage_distribution.svg", format = "svg", height = 10, width = 10, units = 'cm', verbose = False)
o_plot.save(filename = s_path_to_the_output_folder + "/gene_coverage_distribution.png", format = "png", height = 10, width = 10, units = 'cm', dpi = 300, verbose = False)




#Теперь, если нужно, делаю графики для ортогрупп

if s_make_orthogroup_coverage_plots == "true":
	#составлю список с покрытиями ортогрупп и список такой же длины, в котором фразы "Single-copy orthogroups" и "Multicopy orthogroups" — по одной фразе для каждого значения покрытия. Эти два списка имеют одинаковую длину.
	f_infile = open(s_path_to_the_folder_with_calculate_AG_results + "/log.txt", "r")
	l_orthogroup_copy_numbers = [] #Значения могут быть "Single-copy orthogroups" или "Multicopy orthogroups"
	l_coverages = [] #медианные покрытия экзонов генов.

	for s_line in f_infile:
		"""
		Для однокопийных ортогрупп беру информацию о покрытии из строк такого вида:
		For the gene contig_1623:64878-86622 from a single-copy orthogroup 72550at7898, the median coverage is 54.0. It was calculated using 1002 positions.
		
		Для многокопийных ортогрупп беру информацию о покрытии из строк такого вида:
		The mean coverage of genes from a multicopy orthogroup 45170at71240 which contains 2 genes is 28.0. It is larger than 21.8, therefore this orthogroup is considered a true multicopy.
		"""
		
		o_regular_expression_results = re.search(r"^For the gene .+ from a single-copy orthogroup .+ the median coverage is ([\d\.\-\+eE]+)\.\s", s_line)
		if o_regular_expression_results:
			n_coverage = float(o_regular_expression_results.group(1))
			
			l_orthogroup_copy_numbers.append("Single-copy\northogroups\n") #В конце названия ставлю \n, чтобы, как и на графике про гены, в подписи было три строчки. Так графики про гены и про ортогруппы гармоничнее смотрятся рядом.
			
			l_coverages.append(n_coverage)
			
		o_regular_expression_results = re.search(r"^The mean coverage of genes from a multicopy orthogroup .+ which contains \d+ genes is ([\d\.\-\+eE]+)\. ", s_line)
		if o_regular_expression_results:
			n_coverage = float(o_regular_expression_results.group(1))
			
			l_orthogroup_copy_numbers.append("Multicopy\northogroups\n") #В конце названия ставлю \n, чтобы, как и на графике про гены, в подписи было три строчки. Так графики про гены и про ортогруппы гармоничнее смотрятся рядом.
			
			l_coverages.append(n_coverage)

	f_infile.close()

	#print(l_orthogroup_copy_numbers)
	#print(l_coverages)

	#конвертирую список в data frame
	o_data_frame = pandas.DataFrame(data = {"": l_orthogroup_copy_numbers, "Sequencing coverage": l_coverages})

	#print(o_data_frame)

	#Максимальное покрытие по вертикальной оси будет также же, как для генов (его я определял выше, оно записано в переменную n_maximum_coverage_to_draw).

	#цвет #20A647. Подобрал этот цвет вручную, как неплохо подходящий для данной задачи. Сначала попробовал экселевский зелёный (#70AD47), но он подходит хуже.
	#plotnine.scales.scale_y_continuous(expand = (0, 0)) я добавляю, чтобы не было отступа снизу по оси ординат.
	#text = plotnine.element_text(linespacing = 1.1) я добавляю потому, что иначе расстояние между строками текста слишком маленькое.

	#print(n_number_of_bins_in_the_sina_plots)

	if s_distribution_area_type == "equal":
		o_plot = plotnine.ggplot(o_data_frame) + plotnine.aes(x = "", y = "Sequencing coverage") + plotnine.geom_sina(color = "#20A647", scale = "area", size = n_point_size, random_state = n_random_seed, bins = n_number_of_bins_in_the_sina_plots	) + plotnine.scales.scale_x_discrete(limits = ["Single-copy\northogroups\n", "Multicopy\northogroups\n"]) + plotnine.scales.scale_y_continuous(expand = (0, 0), limits = (0, n_maximum_coverage_to_draw)) + plotnine.theme_light() + plotnine.theme(text = plotnine.element_text(linespacing = 1.1))
	if s_distribution_area_type == "proportional":
		o_plot = plotnine.ggplot(o_data_frame) + plotnine.aes(x = "", y = "Sequencing coverage") + plotnine.geom_sina(color = "#20A647", scale = "count", size = n_point_size, random_state = n_random_seed, bins = n_number_of_bins_in_the_sina_plots) + plotnine.scales.scale_x_discrete(limits = ["Single-copy\northogroups\n", "Multicopy\northogroups\n"]) + plotnine.scales.scale_y_continuous(expand = (0, 0), limits = (0, n_maximum_coverage_to_draw)) + plotnine.theme_light() + plotnine.theme(text = plotnine.element_text(linespacing = 1.1))
		
	#Делаю картинки.
	#"verbose = False" использую, чтобы Plotnine не писал warnings про то, что он сохраняет картинки. Почему-то когда Plotnine сохраняет картинки он всегда об этом пишет, используя слово "Warning". Это может смутить пользователя.
	o_plot.save(filename = s_path_to_the_output_folder + "/orthogroup_coverage_distribution.svg", format = "svg", height = 10, width = 10, units = 'cm', verbose = False)
	o_plot.save(filename = s_path_to_the_output_folder + "/orthogroup_coverage_distribution.png", format = "png", height = 10, width = 10, units = 'cm', dpi = 300, verbose = False)





