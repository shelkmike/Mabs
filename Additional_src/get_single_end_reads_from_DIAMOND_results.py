#!/usr/bin/env python3
# coding=utf-8

"""
Этот скрипт берёт результаты выравнивания длинных ридов к белкам DIAMOND и выписывает из файла со всеми ридами только те, которые выровнялись. DIAMOND должен был быть запущен с опцией "--outfmt 6 qseqid qlen sseqid evalue bitscore".

Риды во входном наборе могут быть как в формате FASTA, так и в формате FASTQ. Как незаархивированные, так и заархивированные в gzip.

Риды в выходном файле будут незаархивированными. Название выходного файла даёт пользователь, так что если исходные риды были в формате FASTA, а пользователь указал, что название выходного файла заканчивается на ".fastq", то получится выходной файл в формате FASTA с расширением ".fastq". Так что пользователю нужно учитывать, какое расширение выходного файла давать.

Использование скрипта:
python3 get_single_end_reads_from_DIAMOND_results.py all_reads.fastq diamond_results.txt reads_that_have_matches_to_proteins.fastq
"""

import sys
import os
import re
import gzip

if re.search(r"\.gz$", sys.argv[1], flags = re.IGNORECASE):
	f_all_reads = gzip.open(sys.argv[1], "rt")
else:
	f_all_reads = open(sys.argv[1], "r")

f_DIAMOND_results = open(sys.argv[2], "r")
f_output_reads = open(sys.argv[3], "w")

s_format_of_reads = "" #формат ридов. Значение: "FASTQ" или "FASTA".
if re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", sys.argv[1], flags = re.IGNORECASE):
	s_format_of_reads = "FASTQ"
else:
	s_format_of_reads = "FASTA"

#сначала сделаю список, в котором значения это заголовки всех ридов, которые имели матчи к белкам BUSCO. Важно, что поскольку DIAMOND я запускал с командой "--outfmt 6 qseqid qlen sseqid evalue bitscore", то заголовки тут только до первого пробела.

l_titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character = []
for s_line in f_DIAMOND_results:
	"""
	ERR2173372.5082 250     23397at71240_2  1.67e-18        77.8
	ERR2173372.5096 250     24266at71240_7  1.93e-49        166
	ERR2173372.5096 250     24266at71240_7  3.77e-32        116
	ERR2173372.5110 250     45913at71240_0  4.98e-17        73.6
	ERR2173372.5181 250     28072at71240_6  3.00e-07        45.8
	"""
	
	o_regular_expression_results = re.search(r"^(\S+)", s_line)
	if o_regular_expression_results:
		l_titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character.append(o_regular_expression_results.group(1))

#делаю из списка сет. Насколько я понял (https://stackoverflow.com/questions/7571635/fastest-way-to-check-if-a-value-exists-in-a-list) поиск по сетам быстрее, чем по спискам.
o_set__titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character = set(l_titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character)

#теперь иду по файлу со всеми ридами и выписываю те, у которых часть, которая в заголовке рида до первого пробельного символа, есть в сете o_set__titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character.
n_line_number = 1 #номер строки. считается от единицы.
s_line = f_all_reads.readline()
s_should_I_print_this_read = "no" #"yes", если рид имел матч к белкам BUSCO. "no", если не имел.
while s_line:
	"""
	@ERR2173372.5001 MISEQ:13:000000000-AE4FK:1:1101:13446:2177 length=250
	TCCTCGATTACACTTCTGTCGCCATTAAAATCTTGAAGTCGGGTATCACAGAGGGACTGAAACAGTTCCAACAAGAGGTACTGTGTTAGTTTATATTGCTCTCCTTATTTTGCTTCTCCTTCTCCTGACTCATACATATGTTTCTTCAATAGATTGAGGTTCTTAGTAGCATGAGGCACCCTAATATGGTAATCCTTCTCGGTGCATGTCCCGAGTATGGTTGTCTTGTATATGAGTATATGGAAAATGG
	+ERR2173372.5001 MISEQ:13:000000000-AE4FK:1:1101:13446:2177 length=250
	CCCCCCCFFCFFGGGGGGGGGGGGGHHHHHHHHHHHHHHHGGGHEHHHHHHHHGGGGHHGHHHHHHHGHHHHGHHHHHGHHHHHHHHHHHHHHHHHHGGHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHFHHFGHHHHHHHHHHHGHHHGFHGGHHHHHHHGHHHHHHHHHHHGGDFGHGHHHHHGGGEFHHHDGGGFHFHHGEHHHBGHHBGHHHEFEBFFF
	"""
	#если это строка с заголовком
	if s_format_of_reads == "FASTQ":
		if (n_line_number - 1) % 4 == 0:
			o_regular_expression_results = re.search(r"^\@(\S+)", s_line)
			
			s_read_title = o_regular_expression_results.group(1)

			if s_read_title in o_set__titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character:
				s_should_I_print_this_read = "yes"
			else:
				s_should_I_print_this_read = "no"
	if s_format_of_reads == "FASTA":

		o_regular_expression_results = re.search(r"^>(\S+)", s_line)
		if o_regular_expression_results:	
			s_read_title = o_regular_expression_results.group(1)

			if s_read_title in o_set__titles_of_reads_that_have_matches_to_BUSCO_proteins__before_the_first_whitespace_character:
				s_should_I_print_this_read = "yes"
			else:
				s_should_I_print_this_read = "no"
	
	if s_should_I_print_this_read == "yes":
		f_output_reads.write(s_line)

	s_line = f_all_reads.readline()
	n_line_number += 1
	