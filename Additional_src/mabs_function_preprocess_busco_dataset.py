#!/usr/bin/env python3
# coding=utf-8

"""
Этот модуль содержит одну функцию function_preprocess_busco_dataset, которая делает следующее:
1) Смотрит консервативность последовательностей во всех ортогруппах BUSCO с помощью программы hmmstat из набора программ HMMER. Более высокая консервативность считается у тех последовательностей, у которых более низкая энтропия на позицию.
2) Оставляет самые консервативные ортогруппы — то количество, которое указал пользователь ключом --number_of_BUSCO_orthogroups . Если пользователь указал больше ортогрупп, чем есть в этом наборе BUSCO, то эта функция напишет warning в логи Mabs, и возьмёт все ортогруппы. 
3) Получившуюся базу данных эта функция записывает в папку "BUSCO_dataset_to_use"
"""

import sys
import os
import subprocess
import re

"""
s_path_to_a_local_busco_dataset — путь к разархивированной папке BUSCO.
s_number_of_busco_orthogroups_to_use — сколько ортогрупп BUSCO использовать. Это строка, содержащая или число, или слово "all", если нужно использовать все.
s_path_to_the_output_folder — путь к выходной папке, внутри которой будет создана папка "BUSCO_dataset_to_use".
f_general_logfile — файл с логами. Он был создан ещё скриптом, который вызвал эту функцию.
"""
def function_preprocess_busco_dataset(s_path_to_a_local_busco_dataset, s_number_of_busco_orthogroups_to_use, s_path_to_the_output_folder, f_general_logfile):
	
	s_path_to_the_Additional_folder = os.path.abspath(os.path.dirname(__file__)) #Путь к папке "Additional".
	
	#смотрю в файле dataset.cfg, на основании скольких ортогрупп был сделан этот набор.
	n_full_number_of_busco_orthogroups_in_the_downloaded_dataset = 0
	f_infile = open(s_path_to_a_local_busco_dataset + "/dataset.cfg", "r")
	for s_line in f_infile:
		#number_of_BUSCOs=2326
		if re.search(r"number_of_BUSCOs=(\d+)", s_line):
			o_regular_expression_results = re.search(r"number_of_BUSCOs=(\d+)", s_line)
			n_full_number_of_busco_orthogroups_in_the_downloaded_dataset = int(o_regular_expression_results.group(1))

	if n_full_number_of_busco_orthogroups_in_the_downloaded_dataset == 0:
		f_logs.write("Error. Cannot find the number of orthogroups used in this BUSCO dataset. Stopping.\n")
		sys.exit("Error. Cannot find the number of orthogroups used in this BUSCO dataset. Stopping.\n")

	f_infile.close()

	#если пользователь указал количество ортогрупп в виде слова "all".
	if s_number_of_busco_orthogroups_to_use == "all":
		os.system("cp -r " + s_path_to_a_local_busco_dataset + " " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use")
			
		#конкатенирую все файлы hmm в один 
		os.system("cat " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms/*.hmm >" + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/concatenated_profile_HMMs_of_orthogroups.hmm")
	
	#если пользователь указал количество ортогрупп не в виде слова "all", а в виде целого числа.
	else:
		n_number_of_busco_orthogroups_to_use = int(s_number_of_busco_orthogroups_to_use)
		#если пользователь указал столько ортогрупп, сколько есть в этом наборе BUSCO, то Mabs берёт все ортогруппы.
		if (n_number_of_busco_orthogroups_to_use == n_full_number_of_busco_orthogroups_in_the_downloaded_dataset):
			os.system("cp -r " + s_path_to_a_local_busco_dataset + " " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use")
			
			#конкатенирую все файлы hmm в один 
			os.system("cat " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms/*.hmm >" + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/concatenated_profile_HMMs_of_orthogroups.hmm")
		#если пользователь указал больше ортогрупп, чем есть в этом наборе BUSCO, то Mabs берёт все ортогруппы, и, вдобавок, пишет Warning в основной файл с логами.
		elif (n_number_of_busco_orthogroups_to_use > n_full_number_of_busco_orthogroups_in_the_downloaded_dataset):
			os.system("cp -r " + s_path_to_a_local_busco_dataset + " " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use")
			f_general_logfile.write("Warning: you directed Mabs to use " + str(n_number_of_busco_orthogroups_to_use) + " BUSCO orthogroups. However, the dataset " + s_path_to_a_local_busco_dataset + " contains only " + str(n_full_number_of_busco_orthogroups_in_the_downloaded_dataset) + " orthogroups. Hence, Mabs will use all " + str(n_full_number_of_busco_orthogroups_in_the_downloaded_dataset) + " orthogroups. This is not an error and will not lead to incorrect results of Mabs.")
			
			#конкатенирую все файлы hmm в один 
			os.system("cat " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms/*.hmm >" + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/concatenated_profile_HMMs_of_orthogroups.hmm")
		
		#если пользователь указал меньше ортогрупп, чем есть в этом наборе BUSCO.
		else:
			#сначала сделаю анализ hmmstat для всех HMM, чтобы понять, какие из них наиболее консервативные. Мои тесты на искусственных последовательностях белков показывают, что средняя позиционная относительная энтропия ("mean positional relative entropy", "p relE") действительно имеет более низкие значения для множественных выравниваний, в которых белки меньше отличаются. А вот "relent" почему-то и для выравнивания консервативных белков, и для выравнивания неконсервативных получался 0.59.
			d_orthogroup_name_to_MPRE = {} #словарь, в котором ключ это название ортогруппы, а значение это средняя позиционная относительная энтропия, "mean positional relative entropy", которую я сокращаю до MPRE.
			
			#При анализе, hmmstat будет писать результаты для всех ортогрупп в один файл s_path_to_the_output_folder + "/analysis_of_BUSCO_HMMs.txt" . Потом я его распаршу и получу, собственно, значения MPRE.
			os.system("rm -f " + s_path_to_the_output_folder + "/analysis_of_BUSCO_HMMs.txt")
			
			l_files_in_the_folder_hmms = os.listdir(s_path_to_a_local_busco_dataset + "/hmms") #список файлов в папке hmms. Они имеют названия вида 177652at71240.hmm .
			for s_filename in l_files_in_the_folder_hmms:
				o_regular_expression_results = re.search(r"(.+)\.hmm$", s_filename)
				if o_regular_expression_results:
					os.system(s_path_to_the_Additional_folder + "/HMMER/src/hmmstat " + s_path_to_a_local_busco_dataset + "/hmms/" + s_filename + " >>" + s_path_to_the_output_folder + "/analysis_of_BUSCO_HMMs.txt")
			
			#Теперь иду по файлу analysis_of_BUSCO_HMMs.txt и заполняю словарь d_orthogroup_name_to_MPRE .
			f_infile = open(s_path_to_the_output_folder + "/analysis_of_BUSCO_HMMs.txt", "r")
			for s_line in f_infile:
				"""
				# hmmstat :: display summary statistics for a profile file
				# HMMER 3.2.1 (June 2018); http://hmmer.org/
				# Copyright (C) 2018 Howard Hughes Medical Institute.
				# Freely distributed under the BSD open source license.
				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

				#
				# idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
				# ---- -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
				1      98778at71240         -                  30     0.50    296   0.59   0.62   0.52   0.07
				"""
				o_regular_expression_results = re.search(r"^\d+\s+(\S+).+\s+([\d\.eE\+\-]+)\s+[\d\.eE\+\-]+$", s_line)
				if o_regular_expression_results:
					s_orthogroup_name = o_regular_expression_results.group(1)
					n_MPRE = float(o_regular_expression_results.group(2))
					d_orthogroup_name_to_MPRE[s_orthogroup_name] = n_MPRE
			f_infile.close()	
			
			l_orthogroup_names_sorted_by_increasing_MPRE = sorted(d_orthogroup_name_to_MPRE, key=d_orthogroup_name_to_MPRE.get)
			
			l_names_of_orthogroups_to_take = l_orthogroup_names_sorted_by_increasing_MPRE[0:n_number_of_busco_orthogroups_to_use] #список с названиями ортогрупп, которые я беру.
			os.system("mkdir " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use") #список с теми ортогруппами, которые нужно использовать.
			
			"""
			Теперь я выполняю 6 операций:
			1) Беру белки из файла ancestral, относящиеся к нужным ортогруппам.
			2) Беру относящиеся к нужным ортогруппам строки из файла lengths_cutoff
			3) Беру относящиеся к нужным ортогруппам строки из файла scores_cutoff
			4) Меняет файл dataset.cfg так, чтобы в нём в строке number_of_BUSCOs=2326 стало число n_number_of_busco_orthogroups_to_use. В принципе, Mabs эту строку не использует, но я сделаю это для порядка.
			5) Из папки hmms оставляет только файлы, соответствующие нужным ортогруппам.
			6) Делает файл, в котором марковские модели (файлы .hmm) объединены в один файл.
			
			Тут я беру только те файлы, которые нужны Mabs, лишние не беру.
			"""

			#Операция 1.
			f_infile = open(s_path_to_a_local_busco_dataset + "/ancestral", "r")
			f_outfile = open(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/ancestral", "w")
			
			s_should_I_take_this_protein = "no" #переменная, в которой "no", если белок, на который я смотрю на данный момент, не принадлежит нужной ортогруппе. Если принадлежит, то "yes".
			for s_line in f_infile:
			#>110468at71240
				o_regular_expression_results = re.search(r"^>(.+)", s_line)
				if o_regular_expression_results:					
					s_orthogroup_name = o_regular_expression_results.group(1)
					if s_orthogroup_name in l_names_of_orthogroups_to_take:
						s_should_I_take_this_protein = "yes"
					else:
						s_should_I_take_this_protein = "no"
				
				if s_should_I_take_this_protein == "yes":
					f_outfile.write(s_line)
			f_infile.close()
			f_outfile.close()

			#Операция 2.
			f_infile = open(s_path_to_a_local_busco_dataset + "/lengths_cutoff", "r")
			f_outfile = open(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/lengths_cutoff", "w")
			
			for s_line in f_infile:
			#110468at71240   0       25.5000316255   237
			#84631at71240    0       31.4    314.0
			#139803at71240   0       21.7    217.0

				if re.search(r"^(\S+)", s_line):				
					o_regular_expression_results = re.search(r"(\S+)", s_line)
					s_orthogroup_name = o_regular_expression_results.group(1)
					if s_orthogroup_name in l_names_of_orthogroups_to_take:
						f_outfile.write(s_line)
					
			f_infile.close()
			f_outfile.close()

			#Операция 3.
			f_infile = open(s_path_to_a_local_busco_dataset + "/scores_cutoff", "r")
			f_outfile = open(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/scores_cutoff", "w")
			
			for s_line in f_infile:
			#110468at71240	126.42
			#84631at71240	378.98
			#139803at71240	125.16

				if re.search(r"^(\S+)", s_line):				
					o_regular_expression_results = re.search(r"(\S+)", s_line)
					s_orthogroup_name = o_regular_expression_results.group(1)
					if s_orthogroup_name in l_names_of_orthogroups_to_take:
						f_outfile.write(s_line)
					
			f_infile.close()
			f_outfile.close()

			#Операция 4.
			f_infile = open(s_path_to_a_local_busco_dataset + "/dataset.cfg", "r")
			f_outfile = open(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/dataset.cfg", "w")
			
			for s_line in f_infile:
			#number_of_BUSCOs=2326
				if re.search(r"number_of_BUSCOs=", s_line):				
					f_outfile.write("number_of_BUSCOs=" + str(n_number_of_busco_orthogroups_to_use) + "\n")
				else:
					f_outfile.write(s_line)
					
			f_infile.close()
			f_outfile.close()

			#Операция 5.
			os.system("mkdir " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms")
			l_files_in_the_folder_hmms = os.listdir(s_path_to_a_local_busco_dataset + "/hmms") #список файлов в папке hmms. Они имеют названия вида 177652at71240.hmm .
			for s_filename in l_files_in_the_folder_hmms:
				if re.search(r"(.+)\.hmm$", s_filename):
					o_regular_expression_results = re.search(r"(.+)\.hmm$", s_filename)
					s_orthogroup_name = o_regular_expression_results.group(1)
					if s_orthogroup_name in l_names_of_orthogroups_to_take:
						os.system("cp -r " + s_path_to_a_local_busco_dataset + "/hmms/" + s_filename + " " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms")

			#Операция 6.
			os.system("cat " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/hmms/*.hmm >" + s_path_to_the_output_folder + "/BUSCO_dataset_to_use/concatenated_profile_HMMs_of_orthogroups.hmm")




