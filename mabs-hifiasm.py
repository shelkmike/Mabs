#!/usr/bin/env python3
# coding=utf-8

"""
Mabs-hifiasm, a part of the genome assembly suite Mabs. See https://github.com/shelkmike/Mabs

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
"""

import sys
import os
import re
import datetime
import urllib.request
#import ssl
import math
import shutil
import subprocess
from Additional import mabs_function_preprocess_busco_dataset

if __name__ == '__main__':
	
	s_path_to_the_folder_where_Mabs_lies = os.path.abspath(os.path.dirname( __file__ )) #Путь к папке, где лежит Mabs. Делаю, как написано на https://csatlas.com/python-script-path/
	
	#Сначала проверяю, все ли нужные программы доступны, а также то, что присутствуют папки "Additional" и "Test_datasets". Все проблемы запишу в список l_unavailable_files_and_folders, и потом напечатаю его. Если пользователь допустил ошибки ещё и в командной строке, то напечатаю оба списка проблем (недоступные файлы и ошибки в командной строке) сразу.
	l_unavailable_files_and_folders = []
	if shutil.which("python3") is None:
		l_unavailable_files_and_folders.append("\"python3\" is not in $PATH")
		
	if shutil.which("perl") is None:
		l_unavailable_files_and_folders.append("\"perl\" is not in $PATH")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Bedtools/bedtools"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Bedtools/bedtools has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond has not been found")
	else:
		#Если DIAMOND есть, то проверю, что он версии как минимум 2.0.0
		s_command_results = subprocess.getoutput(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond --version") #Строка имеет вид "diamond version 2.0.13"
		if re.search(r" [01]", s_command_results):
			l_unavailable_files_and_folders.append("The version of DIAMOND should be at least 2.0.0")
			
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmstat"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmstat has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmsearch"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmsearch has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/MetaEuk/metaeuk"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/MetaEuk/metaeuk has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2 has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py"):
		l_unavailable_files_and_folders.append("The file get_single_end_reads_from_DIAMOND_results.py should be in the subfolder \"Additional\" of the folder where Mabs lies.")
	
	if not os.path.isdir(s_path_to_the_folder_where_Mabs_lies + "/Test_datasets"):
		l_unavailable_files_and_folders.append("The subfolder \"Test_datasets\" should be in the folder where Mabs lies.")
	
	
	#делаю парсинг аргументов командной строки. Можно было бы использовать argparse, но когда я делаю это без библиотек, то больше возможностей для того, чтобы сделать интерфейс таким, какой мне нравится.

	s_command_line = " ".join(sys.argv) #команда, которой запущен Mabs-hifiasm, в одну строку.
	s_command_line_reduced = s_command_line #то же, что s_command_line, но после того, как я распаршу какой-нибудь аргумент, я удалю его из этой строки. Если останется какой-то нераспарсенный аргумент, значит пользователь ввёл неизвестные Mabs-hifiasm аргументы, и нужно выдать ошибку.

	#инициализирую исходные значения переменных
	s_path_to_pacbio_hifi_reads = "" #путь к файлу с ридами PacBio HiFi.
	s_path_to_ultralong_nanopore_reads = "" #путь к файлу с ультрадлинными ридами Нанопора.
	s_path_to_hic_short_reads_R1 = "" #путь к файлу с короткими ридами HiC первого конца.
	s_path_to_hic_short_reads_R2 = "" #путь к файлу с короткими ридами HiC второго конца.
	s_busco_dataset_name_online = "" #название файла с базой данных BUSCO с сайта http://mikeshelk.site/Data/BUSCO_datasets/Latest/ (будет непустым только если пользователь использовал опцию "--download_busco_dataset") 
	s_path_to_a_local_busco_dataset = "" #путь к архивированному gzip файлу с датасетом BUSCO на диске или разархивированной папке с датасетом BUSCO на диске.
	n_number_of_cpu_threads_to_use = 10 #количество ядер, которые будет использовать Mabs-hifiasm.
	s_path_to_the_output_folder = "./Mabs_results" #путь к выходной папке Mabs-hifiasm.
	s_genome_size_estimate = "auto" #оценка размера генома.
	n_ploidy = 2 #плоидность.

	s_number_of_busco_orthogroups_to_use = "1000" #сколько ортогрупп BUSCO использовать. Это строка, содержащая или число, или слово "all", если нужно использовать все. Если пользователь укажет больше, чем есть в используемой базе данных BUSCO, то Mabs-hifiasm всё равно будет использовать все.
	s_maximum_allowed_intron_length = "from_BUSCO" #максимальная разрешённая длина интрона. По умолчанию, используется значение из файла dataset.cfg датасета BUSCO.
	
	s_additional_hifiasm_parameters = "" #дополнительные параметры Hifiasm.
	
	s_Mabs_version = "2.18"

	l_errors_in_command_line = [] #список ошибок в командной строке. Если пользователь совершил много ошибок, то Mabs-hifiasm напишет про них все, а не только про первую встреченную.

	#если нет ни одного аргумента командной строки, или есть аргумент командной строки --help, то печатаю хелп
	if (len(sys.argv) == 1) or re.search(r"\s\-\-help", s_command_line):
		print("""Mabs-hifiasm, a program for genome assembly.

Main options:
1) --pacbio_hifi_reads       Path to PacBio HiFi reads, also known as CCS reads.
2) --ultralong_nanopore_reads        Path to ultralong Nanopore reads.
3) --short_hi-c_reads_R1     Path to short Hi-C reads, first reads from a pair.
4) --short_hi-c_reads_R2     Path to short Hi-C reads, second reads from a pair.

[any of the above files may be in FASTQ or FASTA, gzipped or not]

5) --download_busco_dataset        Name of a file from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ . It should be the most taxonomically narrow dataset for your species. For example, for a human genome assembly, use "--download_busco_dataset primates_odb10.2021-02-19.tar.gz" and for a drosophila genome assembly use "--download_busco_dataset diptera_odb10.2020-08-05.tar.gz". Mabs-hifiasm will download the respective file. This option is mutually exclusive with "--local_busco_dataset".
6) --threads        Number of CPU threads to be used by Mabs-hifiasm. The default value is 10.
7) --output_folder        Output folder for Mabs-hifiasm results. The default is "Mabs_results".
8) --number_of_busco_orthogroups        How many BUSCO orthogroups should Mabs-hifiasm use. Should be either a positive integral value or "all" to use all orthogroups. The default value is 1000. 
9) --ploidy		Ploidy of the genome. The default value is 2.
10) --genome_size		Haploid genome size. Should be either "auto" for automatic estimation, or a number, possibly ending with "k", "m" or "g". For example, 1.5g means 1.5 gigabases. The default value is "auto".
11) --max_intron_length        Maximum allowed length of an intron. Should be either "from_BUSCO" to use a value from a BUSCO dataset, or a number, possibly ending with "k", "m" or "g". For example, 20k means 20 kilobases. The default is "from_BUSCO". Change --max_intron_length if you assemble a genome with unusually long introns.
12) --local_busco_dataset        Path to a local BUSCO dataset, manually pre-downloaded from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ or http://busco-data.ezlab.org/v5/data/lineages/. Example: "--local_busco_dataset /home/test/Data/primates_odb10.2021-02-19.tar.gz". May be a .tar.gz file or a decompressed folder. This option is mutually exclusive with "--download_busco_dataset".
13) --additional_hifiasm_parameters        A string with additional parameters to be passed to Hifiasm, enclosed in square brackets. Example: "--additional_hifiasm_parameters [--hom-cov 50 -1 pat.yak -2 mat.yak]".

Informational options:
14) --help        Print this help.
15) --version        Print the version of Mabs-hifiasm.
16) --run_test        Assemble a small dataset to test whether Mabs-hifiasm was installed properly. The assembly takes approximately 10 minutes. If a non-empty file ./Mabs_results/The_best_assembly/assembly.fasta appears, then Mabs-hifiasm was installed successfully.

Example 1:
mabs-hifiasm.py --pacbio_hifi_reads hifi_reads.fastq --download_busco_dataset eudicots_odb10.2020-09-10.tar.gz --threads 40

Example 2:
mabs-hifiasm.py --pacbio_hifi_reads hifi_reads.fastq --short_hi-c_reads_R1 hi-c_reads_trimmed_R1.fastq --short_hi-c_reads_R2 hi-c_reads_trimmed_R2.fastq --ultralong_nanopore_reads nanopore_reads.fastq --download_busco_dataset diptera_odb10.2020-08-05.tar.gz --threads 40
		""")
		sys.exit()

	#если пользователь дал Mabs-hifiasm только параметр --run_test, то делаю сборку по тестовому набору ридов.
	if (len(sys.argv) == 2) and re.search(r"\s\-\-run_test", s_command_line):
		s_path_to_pacbio_hifi_reads = s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/pacbio_hifi_test_reads.fastq.gz"
		s_busco_dataset_name_online = "saccharomycetes_odb10.2020-08-05.tar.gz"
		
		#проверяю, были ли недоступны какие-то программы, которые нужны Mabs.
		if len(l_unavailable_files_and_folders) != 0:
			#Если ошибка была всего одна.
			if len(l_unavailable_files_and_folders) == 1:
				print("There was an error with unavailable files:")
				print(l_unavailable_files_and_folders[0])
				sys.exit()
			#Если было больше одной ошибки.
			if len(l_unavailable_files_and_folders) > 1:
				print("There were errors with unavailable files:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_unavailable_files_and_folders:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_unavailable_files_and_folders[n_error_number - 1])
				sys.exit()
		
		#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка.
		if os.path.exists(s_path_to_the_output_folder):
			if len(os.listdir(s_path_to_the_output_folder)) != 0:
				print("Mabs-hifiasm has stopped because the output folder already exists and is not empty.")
				sys.exit()
		
		#создаю выходную папку Mabs
		os.makedirs(s_path_to_the_output_folder, exist_ok = True)
		
		#Качаю базу данных.
		#отключаю проверку сертификата на случай, если сертификат случайно будет просроченным. Иначе скачивания не будет. Один раз такое было с сайтом BUSCO. Решение взято с https://stackoverflow.com/questions/43204012/how-to-disable-ssl-verification-for-urlretrieve
		#ssl._create_default_https_context = ssl._create_unverified_context
		
		s_path_to_a_local_busco_dataset = s_path_to_the_output_folder + "/" + s_busco_dataset_name_online #путь к месту, где будет лежать скачанный архивированный gzip файл с датасетом BUSCO.
		
		#проверяю, доступен ли адрес http://mikeshelk.site/Data/BUSCO_datasets/Latest/. Он может быть недоступен из-за каких-то проблем с сервером. Если не доступен, то рекомендую пользователю скачать базу с http://busco-data.ezlab.org/v5/data/lineages/ и использовать опцию --local_busco_dataset. Проверку делаю примерно как написано на https://stackoverflow.com/questions/1949318/checking-if-a-website-is-up-via-python . А если доступен, то делаю ещё одну проверку — на то, есть ли нужный файл в папке http://mikeshelk.site/Data/BUSCO_datasets/Latest/
		try:
			s_dummy_variable = urllib.request.urlopen("http://mikeshelk.site/Data/BUSCO_datasets/Latest/").getcode()
			#проверяю, доступен ли нужный файл, и если доступен, то качаю его.
			try:
				urllib.request.urlretrieve("http://mikeshelk.site/Data/BUSCO_datasets/Latest/" + s_busco_dataset_name_online, s_path_to_a_local_busco_dataset)
			except:
				l_errors_in_command_line.append("The file " + s_busco_dataset_name_online + " does not exist at http://mikeshelk.site/Data/BUSCO_datasets/Latest/ .")

		except:
			l_errors_in_command_line.append("Unfortunately, http://mikeshelk.site/Data/BUSCO_datasets/Latest/ is currently not accessible. To test Mabs-hifiasm, download the file http://busco-data.ezlab.org/v5/data/lineages/saccharomycetes_odb10.2020-08-05.tar.gz and run the following command:\nmabs-hifiasm.py --pacbio_hifi_reads [PATH TO THE FOLDER WITH MABS]/Test_datasets/pacbio_hifi_test_reads.fastq.gz --local_busco_dataset saccharomycetes_odb10.2020-08-05.tar.gz")
		
		
		if len(l_errors_in_command_line) != 0:
			#Если ошибка была всего одна.
			if len(l_errors_in_command_line) == 1:
				print("There was an error in the command line of Mabs-hifiasm:")
				print(l_errors_in_command_line[0])
			#Если было больше одной ошибки.
			if len(l_errors_in_command_line) > 1:
				print("There were errors in the command line of Mabs-hifiasm:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_errors_in_command_line:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
			
			#Печатаю пустую строку, как разделитель
			print("")
		
		if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
			print("Mabs-hifiasm has stopped.")
			
			sys.exit()
	
		
	#если пользователь запросил у Mabs-hifiasm не сборку тестового набора.
	else:
		#смотрю, запросил ли пользователь версию Mabs-hifiasm
		if (len(sys.argv) == 1) or re.search(r"\s\-\-version", s_command_line):
			print("Mabs-hifiasm " + s_Mabs_version)
			sys.exit()

		#смотрю, указал ли пользователь в командной строке дополнительные параметры Hifiasm. Этот парсинг я делаю перед парсингом опций собственно Mabs-hifiasm, иначе, если у Hifiasm и Mabs-hifiasm есть одноимённые опции, и пользователь укажет одну из таких опций через "--additional_hifiasm_parameters", то Mabs-hifiasm подумает, что это его опция.
		o_regular_expression_results = re.search(r" --additional_hifiasm_parameters \[(.*?)\]", s_command_line_reduced)
		if o_regular_expression_results:
			s_additional_hifiasm_parameters = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#проверяю, что пользователь не дал опцией --additional_hifiasm_parameters следующие опции: -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2 . Это потому, что Mabs-hifiasm их и так использует.
		if re.search(r"\-s ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"-s\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-o ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"-o\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-t ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"-t\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-n\-hap ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--n-hap\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-hg-size ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--hg-size\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-only-primary ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--only-primary\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-ul ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--ul\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-h1 ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--h1\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")
		if re.search(r"\-\-h2 ", s_additional_hifiasm_parameters):
			l_errors_in_command_line.append("You have given Mabs-hifiasm the option \"--h2\" via the option \"--additional_hifiasm_parameters\". The following options cannot be passed via \"--additional_hifiasm_parameters\": -s, -o, -t, --n-hap, --hg-size, --only-primary, --ul, --h1, --h2.")		
		
		#смотрю, дал ли пользователь риды PacBio HiFi
		o_regular_expression_results = re.search(r" --pacbio_hifi_reads (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_pacbio_hifi_reads = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_pacbio_hifi_reads):
				l_errors_in_command_line.append("The file with PacBio HiFi reads " + s_path_to_pacbio_hifi_reads + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
				
		#смотрю, дал ли пользователь ультрадлинные риды Нанопора
		o_regular_expression_results = re.search(r" --ultralong_nanopore_reads (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_ultralong_nanopore_reads = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_ultralong_nanopore_reads):
				l_errors_in_command_line.append("The file with ultralong Nanopore reads " + s_path_to_ultralong_nanopore_reads + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, дал ли пользователь риды Hi-C первого конца
		o_regular_expression_results = re.search(r" --short_hi-c_reads_R1 (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_hic_short_reads_R1 = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_hic_short_reads_R1):
				l_errors_in_command_line.append("The file with short Hi-C reads R1 " + s_path_to_hic_short_reads_R1 + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, дал ли пользователь риды Hi-C второго конца
		o_regular_expression_results = re.search(r" --short_hi-c_reads_R2 (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_hic_short_reads_R2 = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_hic_short_reads_R2):
				l_errors_in_command_line.append("The file with short Hi-C reads R2 " + s_path_to_hic_short_reads_R2 + " does not exist.")
				
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#проверяю, нет ли такого, что пользователь дал короткие риды первого конца, но не дал короткие риды второго.
		if (s_path_to_hic_short_reads_R1 != "") and (s_path_to_hic_short_reads_R2 == ""):
			l_errors_in_command_line.append("You have provided a file with short Hi-C reads R1 but no file with short Hi-C reads R2.")
		
		#проверяю, нет ли такого, что пользователь дал короткие риды второго конца, но не дал короткие риды первого.
		if (s_path_to_hic_short_reads_R1 == "") and (s_path_to_hic_short_reads_R2 != ""):
			l_errors_in_command_line.append("You have provided a file with short Hi-C reads R2 but no file with short Hi-C reads R1.")
		
		#создаю папку, в которой будут результаты Mabs
		o_regular_expression_results = re.search(r" --output_folder (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_the_output_folder = o_regular_expression_results.group(1)
			
			#Если в начале s_path_to_the_output_folder не стоит "." или "/", то, видимо, пользователь имеет в виду подпапку текущей папки, но не указал "./" в начале. В таком случае я добавлю "./" в начало, иначе могут быть проблемы. Не уверен насчёт mabs-flye и mabs-hifiasm, но у calculate_AG это вызывало проблемы.
			if (not re.search(r"^\.", s_path_to_the_output_folder)) and (not re.search(r"^\/", s_path_to_the_output_folder)):
				s_path_to_the_output_folder = "./" + s_path_to_the_output_folder
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка. Не записываю эту ошибку в список l_errors_in_command_line , а сразу останавливаю работу, потому что если выходная папка уже существует, то в неё нельзя качать файлы BUSCO.
		if os.path.exists(s_path_to_the_output_folder):
			if len(os.listdir(s_path_to_the_output_folder)) != 0:
				print("Mabs-hifiasm has stopped because the output folder already exists and is not empty.")
				sys.exit()
		
		os.makedirs(s_path_to_the_output_folder, exist_ok = True)
		
		#проверяю, что пользователь не указал одновременно --download_busco_dataset и --local_busco_dataset
		if re.search(r" --download_busco_dataset (\S+)", s_command_line_reduced) and re.search(r" --local_busco_dataset (\S+)", s_command_line_reduced):
			l_errors_in_command_line.append("You have provided both \"--download_busco_dataset\" and \"--local_busco_dataset\". Only one of these options should be provided.")
		
		#проверяю, что пользователь указал хотя бы одну из опций --download_busco_dataset и --local_busco_dataset.
		if (not re.search(r" --download_busco_dataset (\S+)", s_command_line_reduced)) and (not re.search(r" --local_busco_dataset (\S+)", s_command_line_reduced)):
			l_errors_in_command_line.append("You should use either \"--download_busco_dataset\" or \"--local_busco_dataset\".")
		
		#если пользователь использовал --download_busco_dataset, то качаю соответствующий файл.
		o_regular_expression_results = re.search(r" --download_busco_dataset (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_busco_dataset_name_online = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
			
			#отключаю проверку сертификата на случай, если сертификат случайно будет просроченным. Иначе скачивания не будет. Один раз такое было с сайтом BUSCO. Решение взято с https://stackoverflow.com/questions/43204012/how-to-disable-ssl-verification-for-urlretrieve 
			#ssl._create_default_https_context = ssl._create_unverified_context
			
			s_path_to_a_local_busco_dataset = s_path_to_the_output_folder + "/" + s_busco_dataset_name_online #путь к месту, где будет лежать скачанный архивированный gzip файл с датасетом BUSCO.
		
			#проверяю, доступен ли адрес http://mikeshelk.site/Data/BUSCO_datasets/Latest/. Он может быть недоступен из-за каких-то проблем с сервером. Если не доступен, то рекомендую пользователю скачать базу с http://busco-data.ezlab.org/v5/data/lineages/ и использовать опцию --local_busco_dataset. Проверку делаю примерно как написано на https://stackoverflow.com/questions/1949318/checking-if-a-website-is-up-via-python . А если доступен, то делаю ещё одну проверку — на то, есть ли нужный файл в папке http://mikeshelk.site/Data/BUSCO_datasets/Latest/
			try:
				s_dummy_variable = urllib.request.urlopen("http://mikeshelk.site/Data/BUSCO_datasets/Latest/").getcode()
				
				#проверяю, доступен ли нужный файл, и если доступен, то качаю его.
				try:
					urllib.request.urlretrieve("http://mikeshelk.site/Data/BUSCO_datasets/Latest/" + s_busco_dataset_name_online, s_path_to_a_local_busco_dataset)
				except:
					l_errors_in_command_line.append("The file " + s_busco_dataset_name_online + " does not exist at http://mikeshelk.site/Data/BUSCO_datasets/Latest/ .")

			except:
				l_errors_in_command_line.append("http://mikeshelk.site/Data/BUSCO_datasets/Latest/ is not accessible. Please, download a BUSCO dataset from http://busco-data.ezlab.org/v5/data/lineages/ and use \"--local_busco_dataset\" instead of \"--download_busco_dataset\".")
		
		#если пользователь использовал --local_busco_dataset
		o_regular_expression_results = re.search(r" --local_busco_dataset (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_a_local_busco_dataset = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
			
			#проверяю, что этот файл или эта папка существует.
			if not os.path.exists(s_path_to_a_local_busco_dataset):
				l_errors_in_command_line.append("The local BUSCO dataset " + s_path_to_a_local_busco_dataset + " does not exist.")
			#проверяю, что если это файл, то его имя заканчивается на .tar.gz
			elif os.path.isfile(s_path_to_a_local_busco_dataset) and (not re.search(r"\.tar\.gz$", s_path_to_a_local_busco_dataset)):
				l_errors_in_command_line.append("The file " + s_path_to_a_local_busco_dataset + " is not a .tar.gz file.")	
		
		#смотрю, указал ли пользователь в командной строке используемое количество потоков.
		o_regular_expression_results = re.search(r" --threads (\d+)", s_command_line_reduced)
		if o_regular_expression_results:
			n_number_of_cpu_threads_to_use = int(o_regular_expression_results.group(1))
							
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, указал ли пользователь в командной строке количество ортогрупп BUSCO, которые нужно использовать.
		o_regular_expression_results = re.search(r" --number_of_busco_orthogroups (\d+|all)", s_command_line_reduced)
		if o_regular_expression_results:
			s_number_of_busco_orthogroups_to_use = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, указал ли пользователь в командной строке плоидность.
		o_regular_expression_results = re.search(r" --ploidy (\d+)", s_command_line_reduced)
		if o_regular_expression_results:
			n_ploidy = int(o_regular_expression_results.group(1))
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, указал ли пользователь в командной строке размер генома. Разрешается три варианта формата: число, число с [kmgKMG] на конце, "auto".
		o_regular_expression_results = re.search(r" --genome_size ([\d\.eE\-\+]+[kmgKMG]?|auto)", s_command_line_reduced)
		if o_regular_expression_results:
			s_genome_size_estimate = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, указал ли пользователь в командной строке максимальную разрешённую длину интронов. Разрешается три варианта формата: число, число с [kmgKMG] на конце, "from_BUSCO".
		o_regular_expression_results = re.search(r" --max_intron_length ([\d\.eE\-\+]+[kmgKMG]?|from_BUSCO)", s_command_line_reduced)
		if o_regular_expression_results:
			s_maximum_allowed_intron_length = o_regular_expression_results.group(1)
							
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
			
			if (re.search(r"^[\d\.eE\-\+]+[kmgKMG]$", s_maximum_allowed_intron_length)):
				o_regular_expression_results = re.search(r"^([\d\.eE\-\+]+)[kmgKMG]$", s_maximum_allowed_intron_length)
				n_the_number = float(o_regular_expression_results.group(1)) #та часть, которую пользователь указал как число. Например для 1.1k это будет 1.1 .
				
				if re.search(r"[kK]$", s_maximum_allowed_intron_length):
					s_maximum_allowed_intron_length = str(int(n_the_number * 1000))
				elif re.search(r"[mM]$", s_maximum_allowed_intron_length):
					s_maximum_allowed_intron_length = str(int(n_the_number * 1000000))
				elif re.search(r"[gG]$", s_maximum_allowed_intron_length):
					s_maximum_allowed_intron_length = str(int(n_the_number * 1000000000))		
		
		#проверяю, что пользователь дал риды HiFi.
		if (s_path_to_pacbio_hifi_reads == ""):
			l_errors_in_command_line.append("Mabs-hifiasm requires very accurate reads, provided via --pacbio_hifi_reads. If you don't have them, use Mabs-flye instead.")
		
		#проверяю, не ввёл ли пользователь какие-то несуществующие опции. Это я определяю по тому, что после того, как я распарсил все команды, в строке s_command_line_reduced осталось что-то, кроме названия исполняемого файла Mabs-hifiasm.
		s_command_line_reduced = re.sub(r"^.*?mabs\-hifiasm(\.py)?\s*", "", s_command_line_reduced)
		if s_command_line_reduced != "":
			l_errors_in_command_line.append("You have provided some options which Mabs-hifiasm doesn't know: " + s_command_line_reduced)
			
		#проверяю, были ли недоступны какие-то программы, которые нужны Mabs, и были ли ошибки в командной строке. Если были какие-то из этих проблем, то пишу об этом и завершаю работу Mabs-hifiasm.
		if len(l_unavailable_files_and_folders) != 0:
			#Если ошибка была всего одна.
			if len(l_unavailable_files_and_folders) == 1:
				print("There was an error with unavailable files:")
				print(l_unavailable_files_and_folders[0])
			#Если было больше одной ошибки.
			if len(l_unavailable_files_and_folders) > 1:
				print("There were errors with unavailable files:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_unavailable_files_and_folders:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_unavailable_files_and_folders[n_error_number - 1])
			
			#Печатаю пустую строку, как разделитель
			print("")
			
		
		if len(l_errors_in_command_line) != 0:
			#Если ошибка была всего одна.
			if len(l_errors_in_command_line) == 1:
				print("There was an error in the command line of Mabs-hifiasm:")
				print(l_errors_in_command_line[0])
			#Если было больше одной ошибки.
			if len(l_errors_in_command_line) > 1:
				print("There were errors in the command line of Mabs-hifiasm:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_errors_in_command_line:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
			
			#Печатаю пустую строку, как разделитель
			print("")
		
		if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
			#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме равно 1
			if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) == 1: 
				print("Mabs-hifiasm has stopped. Please, fix this error and restart Mabs-hifiasm.")
			#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме больше 1
			if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) > 1:
				print("Mabs-hifiasm has stopped. Please, fix these errors and restart Mabs-hifiasm.")
			
			sys.exit()

	f_logs = open(s_path_to_the_output_folder + "/mabs_logs.txt","w",buffering=1) #f_logs это общий файл с логами Mabs-hifiasm, в отличие от трёх дополнительных файлов с логами, которые ведут три отдельных экземпляра Mabs-hifiasm. buffering=1 означает, что буферизация идёт только на уровне строк.
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Started Mabs-hifiasm\n\n")

	f_logs.write("You have run Mabs-hifiasm of version " + s_Mabs_version + " with the following command: " + s_command_line + "\n\n")
	
	#Это строка, в которой указаны пути ко всем ридам, которые нужно давать Modified_hifiasm, а также, если пользователь указал размер генома, то и размер генома. Например, "--hg-size 1g --h1 hic_reads_R1.fastq --hi2 hic_reads_R1.fastq --ul nanopore_reads.fastq hifi_reads.fastq", если был указан размер генома, и были указаны и риды Hi-C, и ультрадлинные риды Нанопора, и риды HiFi. Или, например, просто "hifi_reads.fastq" если размер генома не был указан, и были только риды HiFi. Эта строка нужна, чтобы Mabs-hifiasm было проще передавать аргументы командной строки Modified_hifiasm. Иначе, передача аргументов несколько осложнена, потому что в зависимости от того, какие опции дал Mabs-hifiasm пользователь, программа Modified_hifiasm нужно передавать разное количество аргументов.
	s_command_line_arguments_with_reads_for_Modified_hifiasm = s_path_to_pacbio_hifi_reads
	
	#Если пользователь дал и ультрадлинные риды Нанопора
	if s_path_to_ultralong_nanopore_reads != "":
		s_command_line_arguments_with_reads_for_Modified_hifiasm = "--ul " + s_path_to_ultralong_nanopore_reads + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm
	
	#Если пользователь дал и риды Hi-C
	if s_path_to_hic_short_reads_R1 != "":
		s_command_line_arguments_with_reads_for_Modified_hifiasm = "--h1 " + s_path_to_hic_short_reads_R1 + " --h2 " + s_path_to_hic_short_reads_R2 + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm
	
	#Если пользователь указал размер генома
	if s_genome_size_estimate != "auto":
		s_command_line_arguments_with_reads_for_Modified_hifiasm = "--hg-size " + s_genome_size_estimate + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm
	
	
	#если пользователь делает сборку тестового набора ридов Mabs-hifiasm, то нужно написать подробности этого тестового набора.
	if (len(sys.argv) == 2) and re.search(r"\s\-\-run_test", s_command_line):
		f_logs.write("As a test, Mabs-hifiasm will assemble the first chromosome of Saccharomyces cerevisiae, which is approximately 200 kbp long, using 40x PacBio HiFi reads.\n\n")
		f_logs.write("The command \"mabs-hifiasm.py --run_test\" is equivalent to the command \"mabs-hifiasm.py --pacbio_hifi_reads " + s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/pacbio_hifi_reads__test_set__for_diploid_assembly.fastq.gz --download_busco_dataset saccharomycetes_odb10.2020-08-05.tar.gz\"\n")
		f_logs.write("If after Mabs-hifiasm finishes you see a file ./Mabs_results/The_best_assembly/assembly.fasta which has a size of approximately 200 kilobytes, then the test succeeded.\n\n")
	
	#если пользователь сказал скачать файл с базой BUSCO или сам дал файл (но не папку), то разархивирую файл и меняю значение переменной s_path_to_a_local_busco_dataset с пути к файлу на путь к папке.
	if os.path.isfile(s_path_to_a_local_busco_dataset):
		s_busco_dataset_name = "" #в этой переменной название базы данных. Например, для файла eudicots_odb10.2020-09-10.tar.gz это будет "eudicots_odb10". После разархивирования файла, путь к которому содержится в переменной s_busco_dataset_archive, получается папка, название которой содержится в переменной s_busco_dataset_name.
		s_busco_dataset_name = re.sub(r"^.+\/", r"", s_path_to_a_local_busco_dataset) 
		s_busco_dataset_name = re.sub(r"\..+", r"", s_busco_dataset_name)
		os.system("tar -zxf " + s_path_to_a_local_busco_dataset + " --directory " + s_path_to_the_output_folder)
		
		s_path_to_a_local_busco_dataset = s_path_to_the_output_folder + "/" + s_busco_dataset_name

	#Оставляю из базы BUSCO только нужное количество (s_number_of_busco_orthogroups_to_use) ортогрупп — тех, которые имеют наиболее консервативные последовательности. Если пользователь указал использовать все ортогруппы, то Mabs-hifiasm использует все. Если пользователь указал больше ортогрупп, чем есть в этом наборе BUSCO, то Mabs-hifiasm использует все и пишет Warning в основной файл с логами.
	mabs_function_preprocess_busco_dataset.function_preprocess_busco_dataset(s_path_to_a_local_busco_dataset, s_number_of_busco_orthogroups_to_use, s_path_to_the_output_folder, f_logs)

	#делаю ссылку на файл "ancestral", давая ему расширение .fasta. Затем делаю базу данных DIAMOND.
	#с помощью os.path.abspath() я получают абсолютный путь. Если он относительный, то это может создать проблемы в работоспособности мягкой ссылки.
	os.symlink(os.path.abspath(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/ancestral"), s_path_to_the_output_folder + "/reference_busco_proteins.fasta")
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond makedb --in " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta -d " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta")

	#теперь определяю, какие риды выравниваются к белкам BUSCO. В будущем, для скорости, при определении покрытия в экзонах генов BUSCO я буду использовать только эти длинные риды.
	#согласно мануалу DIAMOND, он может делать выравнивание, когда query в форматах .fasta, .fasta.gz. .fastq, .fastq.gz.
	#поставил --max-target-seqs 1 --max-hsps 1 , потому что если рид дал хоть одно выравнивание, то я считаю, что он, возможно, относится к гену BUSCO.
	#Формат выдачи довольно простой (6 qseqid qlen sseqid evalue bitscore), потому что я исключил все параметры, которые на https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options помечены звёздочкой, чтобы ускорить выравнивание.
	#0.001 будет давать ложные выравнивания только для каждого тысячного рида, что приемлемо (не очень увеличит время картирования minimap2), зато увеличивает чувствительность.


	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond blastx --query " + s_path_to_pacbio_hifi_reads + " --threads " + str(n_number_of_cpu_threads_to_use) + " --frameshift 15 --sensitive --db " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta --out " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_hifi_reads_to_busco_proteins.txt --outfmt 6 qseqid qlen sseqid evalue bitscore --max-target-seqs 1 --max-hsps 1 --min-orf 1 --evalue 0.001 1>" + s_path_to_the_output_folder + "/diamond__pacbio_hifi_reads_to_busco_proteins__stdout.txt 2>" + s_path_to_the_output_folder + "/diamond__pacbio_hifi_reads_to_busco_proteins__stderr.txt")
	
	#Теперь смотрю, какие длинные риды прибластились к белкам BUSCO, и выписываю последовательности этих ридов.
	#Сначала нужно определить, риды в формате FASTQ или в формате FASTA. От этого зависит расширение выходного файла.
	if re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE):
		s_output_extension = "fastq"
	else:
		s_output_extension = "fasta"
	
	os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py " + s_path_to_pacbio_hifi_reads + " " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_hifi_reads_to_busco_proteins.txt " + s_path_to_the_output_folder + "/pacbio_hifi_reads_that_have_matches_to_busco_proteins." + s_output_extension)
	
	s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/pacbio_hifi_reads_that_have_matches_to_busco_proteins." + s_output_extension

	#Теперь, собственно, начинаю проверку 10 точек методом золотого сечения. Параметр, который я оптимизирую, это параметр "-s" Hifiasm. Стартовый интервал -s: [0;1]. n_point_1 это самая левая в данный момент точка (то есть, с наименьшим -s), n_point_4 это самая правая (то есть, с наибольшим -s), а n_point_2 и n_point_3 это две промежуточные, положение которых, собственно, и определяется золотым сечением.
	n_point_1 = 0 #Нижняя граница пробуемых -s.
	n_point_4 = 1 #Верхняя граница пробуемых -s.
	n_point_2 = round(n_point_1 + ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(n_point_4 - n_point_1), 3) #округлю до третьего знака после запятой, иначе у Питона иногда вылезают числа вроде 0.144200000001
	n_point_3 = round(n_point_4 - ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(n_point_4 - n_point_1), 3)

	#Для 0 и 1 (двух крайних точек стартового интервала) я не делаю измерений, потому что метод золотого сечения этого не требует.
	
	#Это список, в который для каждого проверенного -s будет записан AG. Ключ это -s, а значение это AG.
	d_s_to_AG = {} #Например, [0.362] = 762.
	
	#Анализирую вторую точку.
	n_number_of_the_point_under_analysis = 1
	n_s = n_point_2

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-hifiasm started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. -s in this point is " + str(n_s) + "\n")
	
	
	os.mkdir(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s))
	
	#Делаю сборку, после чего конвертирую файл p_ctg.gfa в FASTA, делая файл assembly.fasta .
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm -s " + str(n_s) + " -o " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly --only-primary --n-hap " + str(n_ploidy) + " -t " + str(n_number_of_cpu_threads_to_use) + " " + s_additional_hifiasm_parameters + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm)
	
	#Название выходного файла зависит от того, давал ли пользователь риды Hi-C или нет.
	#если пользователь не дал риды Hi-C
	if (s_path_to_hic_short_reads_R1 == ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.bp.p_ctg.gfa"
	#если пользователь дал риды Hi-C.
	if (s_path_to_hic_short_reads_R1 != ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.hic.p_ctg.gfa"
	
	#теперь из файла GFA с первичными контигами делаю файл FASTA с ними.
	f_infile = open(s_path_to_gfa_with_primary_contigs, "r")
	f_outfile = open(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta", "w")
	for s_line in f_infile:
		#S       ptg000001l      AGTTTACGTTGAACAACCTCCAGGGTTTGT...
		o_regular_expression_results = re.search(r"^[sS]\s+(\S+)\s+(\S+)", s_line)
		if o_regular_expression_results:
			f_outfile.write(">" + o_regular_expression_results.group(1) + "\n" + o_regular_expression_results.group(2) + "\n")
	f_infile.close()
	f_outfile.close()
	
	s_path_to_the_last_assembly_folder = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/" #путь к последней папке со сборкой. Нужен, чтобы из неё перемещать файлы с расширениями .bin и .utg в новую папку со сборкой. Их присутствие ускоряет сборку.
	
	#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
	os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + " --assembly " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta --pacbio_hifi_reads " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe false --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

	#Беру AG, посчитанный скриптом calculate_AG.py
	if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt"):
		f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt", "r")
		s_line_1 = f_infile.readline()
		#AG is 487
		o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
		n_AG_for_point_2 = int(o_regular_expression_results.group(1))
	else:
		f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
		sys.exit()
	
	d_s_to_AG[n_s] = n_AG_for_point_2
	
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("AG for -s " + str(n_s) + " is " + str(n_AG_for_point_2) + "\n\n")

	#Анализирую третью точку.
	n_number_of_the_point_under_analysis += 1
	n_s = n_point_3

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-hifiasm started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. -s in this point is " + str(n_s) + "\n")
	
	os.mkdir(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s))
	
	#Делаю сборку, после чего конвертирую файл p_ctg.gfa в FASTA, делая файл assembly.fasta .
	
	#Перемещаю из прошлой папки со сборкой в эту файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
	os.system("mv " + s_path_to_the_last_assembly_folder + "/*.bin " + s_path_to_the_last_assembly_folder + "/*utg* " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/")
	
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm -s " + str(n_s) + " -o " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly --only-primary --n-hap " + str(n_ploidy) + " -t " + str(n_number_of_cpu_threads_to_use) + " " + s_additional_hifiasm_parameters + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm)
	
	#Название выходного файла зависит от того, давал ли пользователь риды Hi-C или нет.
	#если пользователь не дал риды Hi-C
	if (s_path_to_hic_short_reads_R1 == ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.bp.p_ctg.gfa"
	#если пользователь дал риды Hi-C.
	if (s_path_to_hic_short_reads_R1 != ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.hic.p_ctg.gfa"
	
	#теперь из файла GFA с первичными контигами делаю файл FASTA с ними.
	f_infile = open(s_path_to_gfa_with_primary_contigs, "r")
	f_outfile = open(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta", "w")
	for s_line in f_infile:
		#S       ptg000001l      AGTTTACGTTGAACAACCTCCAGGGTTTGT...
		o_regular_expression_results = re.search(r"^[sS]\s+(\S+)\s+(\S+)", s_line)
		if o_regular_expression_results:
			f_outfile.write(">" + o_regular_expression_results.group(1) + "\n" + o_regular_expression_results.group(2) + "\n")
	f_infile.close()
	f_outfile.close()

	s_path_to_the_last_assembly_folder = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/" #путь к последней папке со сборкой. Нужен, чтобы перемещать из неё в новую папку со сборкой файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
	
	#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
	os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + " --assembly " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta --pacbio_hifi_reads " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe false --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

	#Беру AG, посчитанный скриптом calculate_AG.py
	if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt"):
		f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt", "r")
		s_line_1 = f_infile.readline()
		#AG is 487
		o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
		n_AG_for_point_3 = int(o_regular_expression_results.group(1))
	else:
		f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
		sys.exit()
	
	d_s_to_AG[n_s] = n_AG_for_point_3
	
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("AG for -s " + str(n_s) + " is " + str(n_AG_for_point_3) + "\n\n")

	#теперь последовательно выбираю остальные 8 точек методом золотого сечения и меряю AG для них.
	while n_number_of_the_point_under_analysis < 10: #"<", а не "<=", потому что увеличение номера точки здесь делается в начале цикла.
		n_number_of_the_point_under_analysis += 1
		
		#Смотрю, какая из двух центральных точек (вторая или третья) имеют меньшее значение AG. Если вторая имеет меньшее ли равное третьей, то выкидываю первую точку и сужаю интервал. Если третья имеет меньшее, чем вторая, то выкидываю четвёртую точку и сужаю интервал. При равных значениях выкидывыю правую. Правую выкидываю потому, что по моим впечатлениям оптимальный -s чаще бывает ближе к 0, чем к 1.
		if n_AG_for_point_2 < n_AG_for_point_3:
			n_point_1 = n_point_2
			n_point_2 = n_point_3
			#n_point_4 не меняется
			n_point_3 = round((n_point_4 - ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(n_point_4 - n_point_1)), 3)
			
			n_AG_for_point_1 = n_AG_for_point_2
			n_AG_for_point_2 = n_AG_for_point_3
			#n_AG_for_point_4 не меняется
			n_AG_for_point_3 = -100 #плейсхолдер. Всё равно это значение я сейчас посчитаю.
			
			#Анализирую третью точку.
			n_s = n_point_3

			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("Mabs-hifiasm started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. -s in this point is " + str(n_s) + "\n")
			
			os.mkdir(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s))
			
			#Делаю сборку, после чего конвертирую файл p_ctg.gfa в FASTA, делая файл assembly.fasta .
			
			#Перемещаю из прошлой папки со сборкой в эту файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
			os.system("mv " + s_path_to_the_last_assembly_folder + "/*.bin " + s_path_to_the_last_assembly_folder + "/*utg* " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/")
			
			os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm -s " + str(n_s) + " -o " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly --only-primary --n-hap " + str(n_ploidy) + " -t " + str(n_number_of_cpu_threads_to_use) + " " + s_additional_hifiasm_parameters + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm)
	
			#Название выходного файла зависит от того, давал ли пользователь риды Hi-C или нет.
			#если пользователь не дал риды Hi-C
			if (s_path_to_hic_short_reads_R1 == ""):
				s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.bp.p_ctg.gfa"
			#если пользователь дал риды Hi-C.
			if (s_path_to_hic_short_reads_R1 != ""):
				s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.hic.p_ctg.gfa"
			
			#теперь из файла GFA с первичными контигами делаю файл FASTA с ними.
			f_infile = open(s_path_to_gfa_with_primary_contigs, "r")
			f_outfile = open(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta", "w")
			for s_line in f_infile:
				#S       ptg000001l      AGTTTACGTTGAACAACCTCCAGGGTTTGT...
				o_regular_expression_results = re.search(r"^[sS]\s+(\S+)\s+(\S+)", s_line)
				if o_regular_expression_results:
					f_outfile.write(">" + o_regular_expression_results.group(1) + "\n" + o_regular_expression_results.group(2) + "\n")
			f_infile.close()
			f_outfile.close()
			
			s_path_to_the_last_assembly_folder = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/" #путь к последней папке со сборкой. Нужен, чтобы перемещать из неё в новую папку со сборкой файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
			
			#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
			os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + " --assembly " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta --pacbio_hifi_reads " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe false --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

			#Беру AG, посчитанный скриптом calculate_AG.py
			if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt"):
				f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt", "r")
				s_line_1 = f_infile.readline()
				#AG is 487
				o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
				n_AG_for_point_3 = int(o_regular_expression_results.group(1))
			else:
				f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
				sys.exit()
			
			d_s_to_AG[n_s] = n_AG_for_point_3
			
			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("AG for -s " + str(n_s) + " is " + str(n_AG_for_point_3) + "\n\n")
			
		elif n_AG_for_point_2 >= n_AG_for_point_3:
			#n_point_1 не меняется
			n_point_4 = n_point_3
			n_point_3 = n_point_2
			n_point_2 = round((n_point_1 + ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(n_point_4 - n_point_1)), 3)
			
			#n_AG_for_point_1 не меняется
			n_AG_for_point_4 = n_AG_for_point_3
			n_AG_for_point_3 = n_AG_for_point_2
			n_AG_for_point_2 = -100 #плейсхолдер. Всё равно это значение я сейчас посчитаю.
			
			#Анализирую вторую точку.
			n_s = n_point_2

			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("Mabs-hifiasm started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. -s in this point is " + str(n_s) + "\n")
			
			os.mkdir(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s))
			
			#Делаю сборку, после чего конвертирую файл p_ctg.gfa в FASTA, делая файл assembly.fasta .
			
			#Перемещаю из прошлой папки со сборкой в эту файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
			os.system("mv " + s_path_to_the_last_assembly_folder + "/*.bin " + s_path_to_the_last_assembly_folder + "/*utg* " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/")
			
			os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm -s " + str(n_s) + " -o " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly --only-primary --n-hap " + str(n_ploidy) + " -t " + str(n_number_of_cpu_threads_to_use) + " " + s_additional_hifiasm_parameters + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm)
	
			#Название выходного файла зависит от того, давал ли пользователь риды Hi-C или нет.
			#если пользователь не дал риды Hi-C
			if (s_path_to_hic_short_reads_R1 == ""):
				s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.bp.p_ctg.gfa"
			#если пользователь дал риды Hi-C.
			if (s_path_to_hic_short_reads_R1 != ""):
				s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.hic.p_ctg.gfa"
			
			#теперь из файла GFA с первичными контигами делаю файл FASTA с ними.
			f_infile = open(s_path_to_gfa_with_primary_contigs, "r")
			f_outfile = open(s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta", "w")
			for s_line in f_infile:
				#S       ptg000001l      AGTTTACGTTGAACAACCTCCAGGGTTTGT...
				o_regular_expression_results = re.search(r"^[sS]\s+(\S+)\s+(\S+)", s_line)
				if o_regular_expression_results:
					f_outfile.write(">" + o_regular_expression_results.group(1) + "\n" + o_regular_expression_results.group(2) + "\n")
			f_infile.close()
			f_outfile.close()
			
			#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
			os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + " --assembly " + s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/assembly.fasta --pacbio_hifi_reads " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe false --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))
			
			s_path_to_the_last_assembly_folder = s_path_to_the_output_folder + "/Assembly_for_-s_" + str(n_s) + "/" #путь к последней папке со сборкой. Нужен, чтобы перемещать из неё в новую папку со сборкой файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
			
			#Беру AG, посчитанный скриптом calculate_AG.py
			if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt"):
				f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_-s_" + str(n_s) + "/AG.txt", "r")
				s_line_1 = f_infile.readline()
				#AG is 487
				o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
				n_AG_for_point_2 = int(o_regular_expression_results.group(1))
			else:
				f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
			
			d_s_to_AG[n_s] = n_AG_for_point_2
			
			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("AG for -s " + str(n_s) + " is " + str(n_AG_for_point_2) + "\n\n")

	#После того, как посчитал AG для всех 10 точек, я смотрю, какая из них дала лучший AG. После этого делаю сборку для этого значения "-s", но на этот раз без использования параметра Modified_hifiasm "--only-primary", потому что тут я хочу сделать все файлы, в том числе файлы с фазированной сборкой — может быть, они будут полезны пользователю. Если две точки дают одинаковый AG, то, для определённости, выбираю ту из них, которая имеет меньший -s.
	n_s_that_provides_maximum_AG = -100
	n_maximum_AG = -100
	for n_s in d_s_to_AG:
		if d_s_to_AG[n_s] > n_maximum_AG:
			n_s_that_provides_maximum_AG = n_s
			n_maximum_AG = d_s_to_AG[n_s]
		
		if (d_s_to_AG[n_s] == n_maximum_AG) and (n_s < n_s_that_provides_maximum_AG):
			n_s_that_provides_maximum_AG = n_s
	
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("The optimal -s is " + str(n_s_that_provides_maximum_AG) + ", it provides AG = " + str(n_maximum_AG) + ". Now performing a full assembly for this value of -s.\n\n")
	os.mkdir(s_path_to_the_output_folder + "/The_best_assembly")
				
	#Перемещаю из прошлой папки со сборкой в эту файлы, названия которых имеют форму *.bin или *utg*. Присутствие этих файлов ускоряет сборку.
	os.system("mv " + s_path_to_the_last_assembly_folder + "/*.bin " + s_path_to_the_last_assembly_folder + "/*utg* " + s_path_to_the_output_folder + "/The_best_assembly/")
	
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Modified_hifiasm/modified_hifiasm -s " + str(n_s_that_provides_maximum_AG) + " -o " + s_path_to_the_output_folder + "/The_best_assembly/assembly --n-hap " + str(n_ploidy) + " -t " + str(n_number_of_cpu_threads_to_use) + " " + s_additional_hifiasm_parameters + " " + s_command_line_arguments_with_reads_for_Modified_hifiasm)
	
	#Название выходного файла зависит от того, давал ли пользователь риды Hi-C или нет.
	#если пользователь не дал риды Hi-C
	if (s_path_to_hic_short_reads_R1 == ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/The_best_assembly/assembly.bp.p_ctg.gfa"
	#если пользователь дал риды Hi-C.
	if (s_path_to_hic_short_reads_R1 != ""):
		s_path_to_gfa_with_primary_contigs = s_path_to_the_output_folder + "/The_best_assembly/assembly.hic.p_ctg.gfa"
	
	#теперь из файла GFA с первичными контигами делаю файл FASTA с ними.
	f_infile = open(s_path_to_gfa_with_primary_contigs, "r")
	f_outfile = open(s_path_to_the_output_folder + "/The_best_assembly/assembly.fasta", "w")
	for s_line in f_infile:
		#S       ptg000001l      AGTTTACGTTGAACAACCTCCAGGGTTTGT...
		o_regular_expression_results = re.search(r"^[sS]\s+(\S+)\s+(\S+)", s_line)
		if o_regular_expression_results:
			f_outfile.write(">" + o_regular_expression_results.group(1) + "\n" + o_regular_expression_results.group(2) + "\n")
	f_infile.close()
	f_outfile.close()
		
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-hifiasm finished. The contigs are in the file " + s_path_to_the_output_folder + "/The_best_assembly/assembly.fasta")





