#!/usr/bin/env python3
# coding=utf-8

"""
Этот скрипт считает auN. Это метрика качества сборки, описанная на https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity

auN выдаётся непосредственно на экран.
Для удобства, округляю auN до ближайшего целого.

Пример:
calculate_auN.py some_file.fasta
"""

import sys
import os
import re

s_path_to_the_input_file = sys.argv[1]

l_lengths_of_sequences = [] #список длин последовательностей.

f_infile = open(s_path_to_the_input_file, "r")

s_current_sequence = "" #последовательность, которую скрипт в данный момент рассматривает.

for s_line in f_infile:
	#если это заголовок последовательности, то скрипт записывает длину прошлой последовательности в список l_lengths_of_sequences и обнуляет значение переменной s_current_sequence
	if re.search("^>", s_line):
		if s_current_sequence != "": #если прошлой последовательности не было, то, скорее всего, это потому что я сейчас смотрю на первую последовательность в файле. Тогда прошлую последовательность не надо добавлять в список l_lengths_of_sequences
			n_sequence_length = len(s_current_sequence)
			l_lengths_of_sequences.append(n_sequence_length)
		s_current_sequence = ""
	else:
		s_current_sequence += re.sub(r"\s", r"", s_line) #убираю пробельные символы, включая символы переноса строки.

#Добавляю длину последней последовательности
n_sequence_length = len(s_current_sequence)
l_lengths_of_sequences.append(n_sequence_length)

#Считаю auN.
n_auN = 0
#Числитель и знаменатель в дроби, по которой считается auN.
n_auN_numerator = 0
n_auN_denominator = 0
for n_sequence_length in l_lengths_of_sequences:
	n_auN_numerator += n_sequence_length * n_sequence_length
	n_auN_denominator += n_sequence_length

#Если не было ни одного контига, то считаю, что auN равен 0.
if n_auN_denominator == 0:
	n_auN = 0
else:
	n_auN = n_auN_numerator / n_auN_denominator

print(str(int(round(n_auN, 0))))
	




