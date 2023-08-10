#!/usr/bin/env python3
# coding=utf-8

"""
Скрипт считает N50 последовательностей в формате FASTA. Выдаёт он только одно число - N50. 

Пример:
python3 calculate_N50.py assembly.fasta
"""

import sys
import os
import re

s_path_to_the_input_file = sys.argv[1]

l_lengths_of_sequences = [] #список длин последовательностей.
n_sum_of_lengths_of_all_sequences = 0 #сумма длин всех последовательностей.


f_infile = open(s_path_to_the_input_file, "r")

s_current_sequence = "" #последовательность, которую скрипт в данный момент рассматривает.

for s_line in f_infile:
	#если это заголовок последовательности, то скрипт записывает длину прошлой последовательности в список l_lengths_of_sequences и обнуляет значение переменной s_current_sequence
	if re.search("^>", s_line):
		if s_current_sequence != "": #если прошлой последовательности не было, то, скорее всего, это потому что я сейчас смотрю на первую последовательность в файле. Тогда прошлую последовательность не надо добавлять в список l_lengths_of_sequences
			n_sequence_length = len(s_current_sequence)
			l_lengths_of_sequences.append(n_sequence_length)
			n_sum_of_lengths_of_all_sequences += n_sequence_length
		s_current_sequence = ""
	else:
		s_current_sequence += re.sub(r"\s", r"", s_line) #убираю пробельные символы, включая символы переноса строки.

#Добавляю длину последней последовательности
n_sequence_length = len(s_current_sequence)
l_lengths_of_sequences.append(n_sequence_length)
n_sum_of_lengths_of_all_sequences += n_sequence_length

#сортирую массив длин последовательностей в обратном порядке.
l_lengths_of_sequences_sorted_backwards = sorted(l_lengths_of_sequences, reverse=True)

#считаю N50.
n_current_sum_of_lengths = 0
for n_sequence_length in l_lengths_of_sequences_sorted_backwards:
	n_current_sum_of_lengths += n_sequence_length
	if n_current_sum_of_lengths >= n_sum_of_lengths_of_all_sequences/2:
		print(str(n_sequence_length))
		sys.exit()
	









