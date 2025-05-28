#!/usr/bin/env python3
# coding=utf-8

"""
Скрипту на вход даётся FASTQ-файл (может быть заархивирован gzip), а скрипт считает точность секвенирования каждого рида, и выдаёт в стандартный вывод медианное значение, имеющее вид, например, "93.8%".
Скрипт не просто усредняет Phred-скоры, а сначала считает точность в каждом нуклеотиде, а потом берёт среднюю точность по риду. Phred score у Нанопора учитывает и вероятность инделов тоже (см. https://labs.epi2me.io/quality-scores/), поэтому получившаяся точность учитывает не только ошибочные замены.

Скрипт считает, что у FASTQ формат Phred+33.

Если скрипту на вход даётся несуществующий файл, то скрипт пишет "The input file with reads does not exists".

Пример использования:
calculate_median_read_accuracy.py some_reads.fastq
"""

import sys
import os
import re
import gzip
import statistics

s_path_to_the_input_file = sys.argv[1]

if os.path.exists(s_path_to_the_input_file):

	if re.search(r"(\.gz)$", s_path_to_the_input_file, flags = re.IGNORECASE):
				f_infile = gzip.open(s_path_to_the_input_file, mode = "rt")
	else:
		f_infile = open(s_path_to_the_input_file, "r")

	#Список, элементы которого это точности ридов, выраженные в процентах. 100%, значит рид идеально точный. В списке по одному элементу на каждый рид.
	l_accuracies_of_reads = []

	n_line_number = 0 #Номер строки. Считается от 1.
	for s_line in f_infile:
		n_line_number += 1
		
		if (n_line_number - 4) % 4 == 0:
			s_quality_line = re.sub(r"[\r\n]+$", "", s_line)
			
			n_accuracy_of_the_read = 0 #Точность рида, в процентах. Сначала посчитаю просто как сумму точностей отдельных нуклеотидов, а потом поделю на длину рида.
			
			for s_character in s_quality_line:
				n_Phred_score_corresponding_to_the_character = ord(s_character) - 33
				n_accuracy_of_the_read += 100 * (1 - (10 ** (- n_Phred_score_corresponding_to_the_character / 10)))
			
			n_accuracy_of_the_read = n_accuracy_of_the_read / len(s_quality_line)
			
			l_accuracies_of_reads.append(n_accuracy_of_the_read)

	f_infile.close()

	n_median_accuracy_of_reads = statistics.median(l_accuracies_of_reads)

	print(str(n_median_accuracy_of_reads) + "%")
else:
	sys.exit("The input file with reads does not exists")