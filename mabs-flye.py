#!/usr/bin/env python3
# coding=utf-8

"""
Mabs-flye, a part of the genome assembly suite Mabs. See https://github.com/shelkmike/Mabs

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
			
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye has not been found")
	else:
		#Если Flye есть, проверю, что он версии как минимум 2.9.1
		s_command_results = subprocess.getoutput(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --version") #Строка имеет вид "2.9.1-b1780"
		o_regular_expression_results = re.search(r"^(\d+)\.(\d+)\.(\d+)", s_command_results)
		if o_regular_expression_results:
			n_first_part_of_the_version = int(o_regular_expression_results.group(1))
			n_second_part_of_the_version = int(o_regular_expression_results.group(2))
			n_third_part_of_the_version = int(o_regular_expression_results.group(3))
			
			if (n_first_part_of_the_version < 2) or ((n_first_part_of_the_version == 2) and (n_second_part_of_the_version < 9)) or ((n_first_part_of_the_version == 2) and (n_second_part_of_the_version == 9) and (n_third_part_of_the_version < 1)):
				l_unavailable_files_and_folders.append("The version of Flye should be at least 2.9.1")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmstat"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmstat has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmsearch"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmsearch has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/MetaEuk/metaeuk"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/MetaEuk/metaeuk has not been found")
		
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2 has not been found")

	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-map"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-map has not been found")
	
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-fix"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-fix has not been found")
	
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk has not been found")
	
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py"):
		l_unavailable_files_and_folders.append("The file get_single_end_reads_from_DIAMOND_results.py should be in the subfolder \"Additional\" of the folder where Mabs lies.")
	
	if not os.path.isdir(s_path_to_the_folder_where_Mabs_lies + "/Test_datasets"):
		l_unavailable_files_and_folders.append("The subfolder \"Test_datasets\" should be in the folder where Mabs lies.")


	#делаю парсинг аргументов командной строки. Можно было бы использовать argparse, но когда я делаю это без библиотек, то больше возможностей для того, чтобы сделать интерфейс таким, какой мне нравится.

	s_command_line = " ".join(sys.argv) #команда, которой запущен Mabs-flye, в одну строку.
	s_command_line_reduced = s_command_line #то же, что s_command_line, но после того, как я распаршу какой-нибудь аргумент, я удалю его из этой строки. Если останется какой-то нераспарсенный аргумент, значит пользователь ввёл неизвестные Mabs-flye аргументы, и нужно выдать ошибку.

	#инициализирую исходные значения переменных
	s_path_to_nanopore_reads = "" #путь к файлу с ридами Нанопора.
	s_path_to_pacbio_hifi_reads = "" #путь к файлу с ридами PacBio HiFi.
	s_path_to_pacbio_clr_reads = "" #путь к файлу с ридами PacBio CLR.
	s_busco_dataset_name_online = "" #название файла с базой данных BUSCO с сайта http://mikeshelk.site/Data/BUSCO_datasets/Latest/ (будет непустым только если пользователь использовал опцию "--download_busco_dataset") 
	s_path_to_a_local_busco_dataset = "" #путь к архивированному gzip файлу с датасетом BUSCO на диске или разархивированной папке с датасетом BUSCO на диске.
	n_number_of_cpu_threads_to_use = 10 #количество ядер, которые будет использовать Mabs-flye.
	s_path_to_the_output_folder = "./Mabs_results" #путь к выходной папке Mabs-flye.
	s_genome_size_estimate = "auto" #оценка размера генома.

	s_number_of_busco_orthogroups_to_use = "1000" #сколько ортогрупп BUSCO использовать. Это строка, содержащая или число, или слово "all", если нужно использовать все. Если пользователь укажет больше, чем есть в используемой базе данных BUSCO, то Mabs-flye всё равно будет использовать все.
	s_maximum_allowed_intron_length = "from_BUSCO" #максимальная разрешённая длина интрона. По умолчанию, используется значение из файла dataset.cfg датасета BUSCO.
	s_additional_flye_parameters = "" #дополнительные параметры Flye.
	
	s_Mabs_version = "2.19"

	l_errors_in_command_line = [] #список ошибок в командной строке. Если пользователь совершил много ошибок, то Mabs-flye напишет про них все, а не только про первую встреченную.

	#если нет ни одного аргумента командной строки, или есть аргумент командной строки --help, то печатаю хелп
	if (len(sys.argv) == 1) or re.search(r"\s\-\-help", s_command_line):
		print("""Mabs-flye, a program for genome assembly.

Main options:
1) --nanopore_reads        Path to Nanopore reads.
2) --pacbio_clr_reads        Path to PacBio CLR reads, also known as "old PacBio" reads.
3) --pacbio_hifi_reads        Path to PacBio HiFi reads, also known as CCS reads.

[any of the above files may be in FASTQ or FASTA, gzipped or not]

4) --download_busco_dataset        Name of a file from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ . It should be the most taxonomically narrow dataset for your species. For example, for a human genome assembly, use "--download_busco_dataset primates_odb10.2021-02-19.tar.gz" and for a drosophila genome assembly use "--download_busco_dataset diptera_odb10.2020-08-05.tar.gz". Mabs-flye will download the respective file. This option is mutually exclusive with "--local_busco_dataset".
5) --threads        Number of CPU threads to be used by Mabs-flye. The default value is 10.
6) --output_folder        Output folder for Mabs-flye results. The default is "Mabs_results".
7) --number_of_busco_orthogroups        How many BUSCO orthogroups should Mabs-flye use. Should be either a positive integral value or "all" to use all orthogroups. The default value is 1000. 
8) --genome_size		Haploid genome size. Should be either "auto" for automatic estimation, or a number ending with "k", "m" or "g". For example, 1.5g means 1.5 gigabases. The default value is "auto".
9) --max_intron_length        Maximum allowed length of an intron. Should be either "from_BUSCO" to use a value from a BUSCO dataset, or a number, possibly ending with "k", "m" or "g". For example, 20k means 20 kilobases. The default is "from_BUSCO". Change --max_intron_length if you assemble a genome with unusually long introns.
10) --local_busco_dataset        Path to a local BUSCO dataset, manually pre-downloaded from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ or http://busco-data.ezlab.org/v5/data/lineages/. Example: "--local_busco_dataset /home/test/Data/primates_odb10.2021-02-19.tar.gz". May be a .tar.gz file or a decompressed folder. This option is mutually exclusive with "--download_busco_dataset".
11) --additional_flye_parameters        A string with additional parameters to be passed to Flye, enclosed in square brackets. Example: "--additional_flye_parameters [--scaffold --min-overlap 20000]".

Informational options:
12) --help        Print this help.
13) --version        Print the version of Mabs-flye.
14) --run_test        Assemble a small dataset to test whether Mabs-flye was installed properly. The assembly takes approximately 10 minutes. If a non-empty file ./Mabs_results/The_best_assembly/assembly.fasta appears, then Mabs-flye was installed successfully.

Example 1:
mabs-flye.py --nanopore_reads nanopore_reads.fastq --download_busco_dataset eudicots_odb10.2020-09-10.tar.gz --threads 40

Example 2:
mabs-flye.py --nanopore_reads nanopore_reads.fastq --pacbio_hifi_reads pacbio_hifi_reads.fastq --download_busco_dataset diptera_odb10.2020-08-05.tar.gz --threads 40
		""")
		sys.exit()

	#если пользователь дал Mabs-flye только параметр --run_test, то делаю сборку по тестовому набору ридов.
	if (len(sys.argv) == 2) and re.search(r"\s\-\-run_test", s_command_line):
		s_path_to_nanopore_reads = s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/nanopore_test_reads.fastq.gz"
		s_path_to_pacbio_clr_reads = s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/pacbio_clr_test_reads.fastq.gz"
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
				print("Mabs-flye has stopped because the output folder already exists and is not empty.")
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
			l_errors_in_command_line.append("Unfortunately, http://mikeshelk.site/Data/BUSCO_datasets/Latest/ is currently not accessible. To test Mabs-flye, download the file http://busco-data.ezlab.org/v5/data/lineages/saccharomycetes_odb10.2020-08-05.tar.gz and run the following command:\nmabs-flye.py --nanopore_reads [PATH TO THE FOLDER WITH MABS]/Test_datasets/nanopore_test_reads.fastq.gz --pacbio_clr_reads [PATH TO THE FOLDER WITH MABS]/Test_datasets/pacbio_clr_test_reads.fastq.gz --local_busco_dataset saccharomycetes_odb10.2020-08-05.tar.gz")
		
		if len(l_errors_in_command_line) != 0:
			#Если ошибка была всего одна.
			if len(l_errors_in_command_line) == 1:
				print("There was an error in the command line of Mabs-flye:")
				print(l_errors_in_command_line[0])
			#Если было больше одной ошибки.
			if len(l_errors_in_command_line) > 1:
				print("There were errors in the command line of Mabs-flye:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_errors_in_command_line:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
			
			#Печатаю пустую строку, как разделитель
			print("")
			
		if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
			print("Mabs-flye has stopped.")
			
			sys.exit()

	#если пользователь запросил у Mabs-flye не сборку тестового набора.
	else:
		#смотрю, запросил ли пользователь версию Mabs-flye
		if (len(sys.argv) == 1) or re.search(r"\s\-\-version", s_command_line):
			print("Mabs-flye " + s_Mabs_version)
			sys.exit()
		
		#смотрю, указал ли пользователь в командной строке дополнительные параметры Flye. Этот парсинг я делаю перед парсингом опций собственно Mabs-flye, иначе, если у Flye и Mabs-flye есть одноимённые опции, и пользователь укажет одну из таких опций через "--additional_hifiasm_parameters", то Mabs-flye подумает, что это его опция.
		o_regular_expression_results = re.search(r" --additional_flye_parameters \[(.*?)\]", s_command_line_reduced)
		if o_regular_expression_results:
			s_additional_flye_parameters = o_regular_expression_results.group(1)
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#проверяю, что пользователь не дал опцией --additional_flye_parameters следующие опции: --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative . Это потому, что Mabs-flye их и так использует.
		if re.search(r"\-\-nano\-raw ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"--nano-raw\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-\-out\-dir ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"--out-dir\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-\-threads ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"--threads\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-\-no\-alt\-contigs ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"--no-alt-contigs\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-\-nano\-raw ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"--genome-size\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-o ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"-o\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-t ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"-t\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"\-g ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"-g\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"assemble_ovlp_divergence ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"assemble_ovlp_divergence\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"repeat_graph_ovlp_divergence ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"repeat_graph_ovlp_divergence\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		if re.search(r"assemble_divergence_relative ", s_additional_flye_parameters):
			l_errors_in_command_line.append("You have given Mabs-flye the option \"assemble_divergence_relative\" via the option \"--additional_flye_parameters\". The following options cannot be passed via \"--additional_flye_parameters\": --nano-raw, --out-dir, --threads, --no-alt-contigs, --genome-size, -o, -t, -g, assemble_ovlp_divergence, repeat_graph_ovlp_divergence, assemble_divergence_relative.")
		
		#смотрю, дал ли пользователь риды Нанопора
		o_regular_expression_results = re.search(r" --nanopore_reads (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_nanopore_reads = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_nanopore_reads):
				l_errors_in_command_line.append("The file with Nanopore reads " + s_path_to_nanopore_reads + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, дал ли пользователь риды PacBio HiFi
		o_regular_expression_results = re.search(r" --pacbio_hifi_reads (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_pacbio_hifi_reads = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_pacbio_hifi_reads):
				l_errors_in_command_line.append("The file with PacBio HiFi reads " + s_path_to_pacbio_hifi_reads + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
		#смотрю, дал ли пользователь риды PacBio CLR
		o_regular_expression_results = re.search(r" --pacbio_clr_reads (\S+)", s_command_line_reduced)
		if o_regular_expression_results:
			s_path_to_pacbio_clr_reads = o_regular_expression_results.group(1)
			if not os.path.isfile(s_path_to_pacbio_clr_reads):
				l_errors_in_command_line.append("The file with PacBio CLR reads " + s_path_to_pacbio_clr_reads + " does not exist.")
			
			s_string_to_remove = re.escape(o_regular_expression_results.group(0))
			s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
		
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
				print("Mabs-flye has stopped because the output folder already exists and is not empty.")
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
		
		#проверяю, что пользователь дал хоть какие-то длинные риды.
		if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads == ""):
			l_errors_in_command_line.append("You have not given Mabs-flye Nanopore reads, PacBio HiFi reads or PacBio CLR reads. Mabs-flye requires at least one set of long reads.")
		
		#проверяю, не ввёл ли пользователь какие-то несуществующие опции. Это я определяю по тому, что после того, как я распарсил все команды, в строке s_command_line_reduced осталось что-то, кроме названия исполняемого файла Mabs-flye.
		s_command_line_reduced = re.sub(r"^.*?mabs\-flye(\.py)?\s*", "", s_command_line_reduced)
		if s_command_line_reduced != "":
			l_errors_in_command_line.append("You have provided some options which Mabs-flye doesn't know: " + s_command_line_reduced)
			
		#проверяю, были ли недоступны какие-то программы, которые нужны Mabs, и были ли ошибки в командной строке. Если были какие-то из этих проблем, то пишу об этом и завершаю работу Mabs-flye.
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
				print("There was an error in the command line of Mabs-flye:")
				print(l_errors_in_command_line[0])
			#Если было больше одной ошибки.
			if len(l_errors_in_command_line) > 1:
				print("There were errors in the command line of Mabs-flye:")
				n_error_number = 0 #порядковый номер ошибки. Считается от 1.
				for s_error_text in l_errors_in_command_line:
					n_error_number += 1
					print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
			
			#Печатаю пустую строку, как разделитель
			print("")
			
		if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
			#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме равно 1
			if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) == 1: 
				print("Mabs-flye has stopped. Please, fix this error and restart Mabs-flye.")
			#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме больше 1
			if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) > 1:
				print("Mabs-flye has stopped. Please, fix these errors and restart Mabs-flye.")
			
			sys.exit()

	f_logs = open(s_path_to_the_output_folder + "/mabs_logs.txt","w",buffering=1) #f_logs это общий файл с логами Mabs-flye, в отличие от трёх дополнительных файлов с логами, которые ведут три отдельных экземпляра Mabs-flye. buffering=1 означает, что буферизация идёт только на уровне строк.
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Started Mabs-flye\n\n")

	f_logs.write("You have run Mabs-flye of version " + s_Mabs_version + " with the following command: " + s_command_line + "\n\n")

	#если пользователь делает сборку тестового набора ридов Mabs-flye, то нужно написать подробности этого тестового набора.
	if (len(sys.argv) == 2) and re.search(r"\s\-\-run_test", s_command_line):
		f_logs.write("As a test, Mabs-flye will assemble the first chromosome of Saccharomyces cerevisiae, which is approximately 200 kbp long, using 20x Nanopore reads and 10x PacBio CLR reads.\n\n")
		f_logs.write("The command \"mabs-flye.py --run_test\" is equivalent to the command \"mabs-flye.py --nanopore_reads " + s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/nanopore_test_reads.fastq.gz --pacbio_clr_reads " + s_path_to_the_folder_where_Mabs_lies + "/Test_datasets/pacbio_clr_test_reads.fastq.gz --download_busco_dataset saccharomycetes_odb10.2020-08-05.tar.gz\"\n")
		f_logs.write("If after Mabs-flye finishes you see a file ./Mabs_results/The_best_assembly/assembly.fasta which has a size of approximately 200 kilobytes, then the test succeeded.\n\n")
	
	#если пользователь сказал скачать файл с базой BUSCO или сам дал файл (но не папку), то разархивирую файл и меняю значение переменной s_path_to_a_local_busco_dataset с пути к файлу на путь к папке.
	if os.path.isfile(s_path_to_a_local_busco_dataset):
		s_busco_dataset_name = "" #в этой переменной название базы данных. Например, для файла eudicots_odb10.2020-09-10.tar.gz это будет "eudicots_odb10". После разархивирования файла, путь к которому содержится в переменной s_busco_dataset_archive, получается папка, название которой содержится в переменной s_busco_dataset_name.
		s_busco_dataset_name = re.sub(r"^.+\/", r"", s_path_to_a_local_busco_dataset) 
		s_busco_dataset_name = re.sub(r"\..+", r"", s_busco_dataset_name)
		os.system("tar -zxf " + s_path_to_a_local_busco_dataset + " --directory " + s_path_to_the_output_folder)
		
		s_path_to_a_local_busco_dataset = s_path_to_the_output_folder + "/" + s_busco_dataset_name
		
			
	#Оставляю из базы BUSCO только нужное количество (s_number_of_busco_orthogroups_to_use) ортогрупп — тех, которые имеют наиболее консервативные последовательности. Если пользователь указал использовать все ортогруппы, то Mabs-flye использует все. Если пользователь указал больше ортогрупп, чем есть в этом наборе BUSCO, то Mabs-flye использует все и пишет Warning в основной файл с логами.
	mabs_function_preprocess_busco_dataset.function_preprocess_busco_dataset(s_path_to_a_local_busco_dataset, s_number_of_busco_orthogroups_to_use, s_path_to_the_output_folder, f_logs)

	#делаю ссылку на файл "ancestral", давая ему расширение .fasta. Затем делаю базу данных DIAMOND.
	#с помощью os.path.abspath() я получают абсолютный путь. Если он относительный, то это может создать проблемы в работоспособности мягкой ссылки.
	os.symlink(os.path.abspath(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/ancestral"), s_path_to_the_output_folder + "/reference_busco_proteins.fasta")
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond makedb --in " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta -d " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta")

	#теперь определяю, какие риды выравниваются к белкам BUSCO. В будущем, для скорости, при определении покрытия на каждой итерации я буду использовать только эти длинные риды.
	#согласно мануалу DIAMOND, он может делать выравнивание, когда query в форматах .fasta, .fasta.gz. .fastq, .fastq.gz.
	#поставил --max-target-seqs 1 --max-hsps 1 , потому что если рид дал хоть одно выравнивание, то я считаю, что он, возможно, относится к гену BUSCO.
	#Формат выдачи довольно простой (6 qseqid qlen sseqid evalue bitscore), потому что я исключил все параметры, которые на https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options помечены звёздочкой, чтобы ускорить выравнивание.
	#0.001 будет давать ложные выравнивания только для каждого тысячного рида, что приемлемо (не очень увеличит время картирования minimap2), зато увеличивает чувствительность.

	#Пути к ридам, которые имеют матчи к белкам BUSCO. Эти переменные (для тех наборов ридов, которые дал пользователь) я заполню ниже.
	s_path_to_nanopore_reads_that_correspond_to_busco_genes = ""
	s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes = ""
	s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes = ""

	#Риды Нанопора.
	if s_path_to_nanopore_reads != "":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond blastx --query " + s_path_to_nanopore_reads + " --threads " + str(n_number_of_cpu_threads_to_use) + " --frameshift 15 --ultra-sensitive --db " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta --out " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_nanopore_reads_to_busco_proteins.txt --outfmt 6 qseqid qlen sseqid evalue bitscore --max-target-seqs 1 --max-hsps 1 --min-orf 1 --evalue 0.001 1>" + s_path_to_the_output_folder + "/diamond__nanopore_reads_to_busco_proteins__stdout.txt 2>" + s_path_to_the_output_folder + "/diamond__nanopore_reads_to_busco_proteins__stderr.txt")
		
		#Теперь смотрю, какие длинные риды прибластились к белкам BUSCO, и выписываю последовательности этих ридов.
		#Сначала нужно определить, риды в формате FASTQ или в формате FASTA. От этого зависит расширение выходного файла.
		if re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_nanopore_reads, flags = re.IGNORECASE):
			s_output_extension = "fastq"
		else:
			s_output_extension = "fasta"
		
		os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py " + s_path_to_nanopore_reads + " " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_nanopore_reads_to_busco_proteins.txt " + s_path_to_the_output_folder + "/nanopore_reads_that_have_matches_to_busco_proteins." + s_output_extension)
		
		s_path_to_nanopore_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/nanopore_reads_that_have_matches_to_busco_proteins." + s_output_extension

	#Риды PacBio HiFi. Для них я использую DIAMOND в режиме не ultra-sensitive, а sensitive.
	if s_path_to_pacbio_hifi_reads != "":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond blastx --query " + s_path_to_pacbio_hifi_reads + " --threads " + str(n_number_of_cpu_threads_to_use) + " --frameshift 15 --sensitive --db " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta --out " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_hifi_reads_to_busco_proteins.txt --outfmt 6 qseqid qlen sseqid evalue bitscore --max-target-seqs 1 --max-hsps 1 --min-orf 1 --evalue 0.001 1>" + s_path_to_the_output_folder + "/diamond__pacbio_hifi_reads_to_busco_proteins__stdout.txt 2>" + s_path_to_the_output_folder + "/diamond__pacbio_hifi_reads_to_busco_proteins__stderr.txt")
		
		#Теперь смотрю, какие длинные риды прибластились к белкам BUSCO, и выписываю последовательности этих ридов.
		#Сначала нужно определить, риды в формате FASTQ или в формате FASTA. От этого зависит расширение выходного файла.
		if re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE):
			s_output_extension = "fastq"
		else:
			s_output_extension = "fasta"
		
		os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py " + s_path_to_pacbio_hifi_reads + " " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_hifi_reads_to_busco_proteins.txt " + s_path_to_the_output_folder + "/pacbio_hifi_reads_that_have_matches_to_busco_proteins." + s_output_extension)
		
		s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/pacbio_hifi_reads_that_have_matches_to_busco_proteins." + s_output_extension
				
	#Риды PacBio CLR.
	if s_path_to_pacbio_clr_reads != "":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond blastx --query " + s_path_to_pacbio_clr_reads + " --threads " + str(n_number_of_cpu_threads_to_use) + " --frameshift 15 --ultra-sensitive --db " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta --out " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_clr_reads_to_busco_proteins.txt --outfmt 6 qseqid qlen sseqid evalue bitscore --max-target-seqs 1 --max-hsps 1 --min-orf 1 --evalue 0.001 1>" + s_path_to_the_output_folder + "/diamond__pacbio_clr_reads_to_busco_proteins__stdout.txt 2>" + s_path_to_the_output_folder + "/diamond__pacbio_clr_reads_to_busco_proteins__stderr.txt")
		
		#Теперь смотрю, какие длинные риды прибластились к белкам BUSCO, и выписываю последовательности этих ридов.
		#Сначала нужно определить, риды в формате FASTQ или в формате FASTA. От этого зависит расширение выходного файла.
		if re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE):
			s_output_extension = "fastq"
		else:
			s_output_extension = "fasta"
		
		os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/Additional/get_single_end_reads_from_DIAMOND_results.py " + s_path_to_pacbio_clr_reads + " " + s_path_to_the_output_folder + "/diamond_results_for_alignment_of_pacbio_clr_reads_to_busco_proteins.txt " + s_path_to_the_output_folder + "/pacbio_clr_reads_that_have_matches_to_busco_proteins." + s_output_extension)
		
		s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/pacbio_clr_reads_that_have_matches_to_busco_proteins." + s_output_extension

	#конкатенирую все длинные риды, относящиеся к генам BUSCO, в один файл. Если все файлы были в формате FASTQ, то и выходной будет FASTQ. Если хотя бы один из входных был в формате FASTA, то и выходной будет в формате FASTA.
	s_path_to_all_long_reads_that_correspond_to_busco_genes = "" #в этой переменной будет путь к файлу со всеми длинными ридами, которые относятся к генам BUSCO.
	#если пользователь дал три набора длинных ридов сразу.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads != ""):
		
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_nanopore_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			os.system("cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + ") > " + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta"
		#если все наборы ридов были в формате FASTQ.
		else:
			os.system("cat " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + " " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + " >" + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq"
		
	#если пользователь дал риды Нанопора и PacBio HiFi.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads == ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_nanopore_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)):
			os.system("cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + ") > " + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta"
		#если все наборы ридов были в формате FASTQ.
		else:
			os.system("cat " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + " " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " >" + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq"

	#если пользователь дал риды Нанопора и PacBio CLR.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads != ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_nanopore_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			os.system("cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + " <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + ") > " + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta"
		#если все наборы ридов были в формате FASTQ.
		else:
			os.system("cat " + s_path_to_nanopore_reads_that_correspond_to_busco_genes + " " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + " >" + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq"

	#если пользователь дал риды PacBio HiFi и PacBio CLR.
	if (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads != ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			os.system("cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + ") > " + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fasta"
		#если все наборы ридов были в формате FASTQ.
		else:
			os.system("cat " + s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes + " " + s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes + " >" + s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq")
			
			s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_the_output_folder + "/all_long_reads_that_correspond_to_busco_genes.fastq"

	#если из длинных ридов пользователь дал только риды Нанопора
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads == ""):
		s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_nanopore_reads_that_correspond_to_busco_genes

	#если из длинных ридов пользователь дал только риды PacBio HiFi
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads == ""):
		s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_pacbio_hifi_reads_that_correspond_to_busco_genes

	#если из длинных ридов пользователь дал только риды PacBio CLR
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads != ""):
		s_path_to_all_long_reads_that_correspond_to_busco_genes = s_path_to_pacbio_clr_reads_that_correspond_to_busco_genes


	#если пользователь дал больше одного набора длинных ридов, то конкатенирую все длинные риды, конвертируя их в формат .fasta.gz. Если все файлы с длинными ридами в формате FASTQ, то Mabs-flye делает конкатенированный файл в формате FASTQ. Если хотя бы один в формате FASTA, то делает конкатенированный файл в формате FASTA.
	#Конкатенацию делаю как написано на https://stackoverflow.com/questions/16013616/os-system-complains-about-round-brackets , потому что если делать просто через os.system, то возникает ошибка, из-за того, что Python пытается использовать sh вместо bash.
	#Устанавливаю значение переменной s_long_reads_option_for_calculate_AG, в которой будет опция, которой эти риды нужно дать calculate_AG.py. Значение может быть "--pacbio_hifi_reads", "--pacbio_clr_reads", "--nanopore_reads". При конкатенации по нескольким типам ридов, я даю их скрипту calculate_AG.py как наименее точные из использованных. Считаю, что HiFi точнее, чем CLR, которые, в свою очередь, точнее, чем риды Нанопора. 
	#если пользователь дал три набора длинных ридов сразу.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads != ""):
		
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fasta.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fasta.gz"
		#если все наборы ридов были в формате FASTQ.
		else:
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_hifi_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fastq.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fastq.gz"
		
		s_long_reads_option_for_calculate_AG = "--nanopore_reads"
		
	#если пользователь дал риды Нанопора и PacBio HiFi.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads == ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_nanopore_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)):
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fasta.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fasta.gz"
		#если все наборы ридов были в формате FASTQ.
		else:
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_hifi_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fastq.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fastq.gz"
		
		s_long_reads_option_for_calculate_AG = "--nanopore_reads"
		
	#если пользователь дал риды Нанопора и PacBio CLR.
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads != ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fasta.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fasta.gz"
		#если все наборы ридов были в формате FASTQ.
		else:
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_nanopore_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fastq.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fastq.gz"
		
		s_long_reads_option_for_calculate_AG = "--nanopore_reads"
		
	#если пользователь дал риды PacBio HiFi и PacBio CLR.
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads != ""):
		#если хотя бы один из наборов ридов был в формате FASTA.
		if (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_hifi_reads, flags = re.IGNORECASE)) or (not re.search(r"(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$", s_path_to_pacbio_clr_reads, flags = re.IGNORECASE)):
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_hifi_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq -a " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fasta.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fasta.gz"
		#если все наборы ридов были в формате FASTQ.
		else:
			subprocess.call(["bash", "-c", "cat <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_hifi_reads + ") <(" + s_path_to_the_folder_where_Mabs_lies + "/Additional/SeqTk/seqtk seq " + s_path_to_pacbio_clr_reads + ") | gzip -1 > " + s_path_to_the_output_folder + "/all_long_reads.fastq.gz"])
			s_path_to_the_file_with_all_long_reads = s_path_to_the_output_folder + "/all_long_reads.fastq.gz"
		
		s_long_reads_option_for_calculate_AG = "--pacbio_clr_reads"
		
	#если из длинных ридов пользователь дал только риды Нанопора
	if (s_path_to_nanopore_reads != "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads == ""):			
		s_path_to_the_file_with_all_long_reads = s_path_to_nanopore_reads
		s_long_reads_option_for_calculate_AG = "--nanopore_reads"
	#если из длинных ридов пользователь дал только риды PacBio HiFi
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads != "") and (s_path_to_pacbio_clr_reads == ""):
		s_path_to_the_file_with_all_long_reads = s_path_to_pacbio_hifi_reads
		s_long_reads_option_for_calculate_AG = "--pacbio_hifi_reads"
	#если из длинных ридов пользователь дал только риды PacBio CLR
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads != ""):					
		s_path_to_the_file_with_all_long_reads = s_path_to_pacbio_clr_reads
		s_long_reads_option_for_calculate_AG = "--pacbio_clr_reads"


	#Теперь, собственно, начинаю проверку 10 точек методом золотого сечения. n_point_1 это самая левая в данный момент точка (то есть, с наименьшим log10(max_divergence)), n_point_4 это самая правая (то есть, с наибольшим log10(max_divergence)), а n_point_2 и n_point_3 это две промежуточные, положение которых, собственно, и определяется золотым сечением.
	#Любые риды я при сборке даю Flye как "nano-raw", то есть нескорректированные риды Нанопора, потому что если мне нужно коллапсировать области генома с очень высокой гетерозиготностью, то для Mabs-flye это примерно эквивалентно тому, что в ридах много ошибок секвенирования. Нужно будет подумать, насколько это правильно. Скрипту calculate_AG.py я тоже даю любые риды как риды Нанопора, то есть с опцией --nanopore_reads.
	n_point_1 = 0.0001 #Нижняя граница пробуемых max_divergence. 0.0001 это 0.01%.
	n_point_4 = 0.5 #Верхняя граница пробуемых max_divergence. 0.5 это 50%.
	n_point_2 = round(10**(math.log10(n_point_1) + ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(math.log10(n_point_4) - math.log10(n_point_1))), 6) #округлю до шестого знака после запятой, иначе у Питона иногда вылезают числа вроде 0.0001442000001
	n_point_3 = round(10**(math.log10(n_point_4) - ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(math.log10(n_point_4) - math.log10(n_point_1))), 6)

	#Для 0.0001 и 0.5 я не делаю измерений, потому что метод золотого сечения этого не требует.
	
	#Это список, в который для каждого проверенного max_divergence будет записан AG. Ключ это max_divergence, а значение это AG.
	d_max_divergence_to_AG = {} #Например, [0.123] = 762.

	#Анализирую вторую точку.
	n_number_of_the_point_under_analysis = 1
	n_max_divergence = n_point_2

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-flye started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. Max_divergence in this point is " + str(n_max_divergence) + "\n")

	#если пользователь не указывал размер генома
	if s_genome_size_estimate == "auto":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
	#если пользователь указал размер генома
	else:
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --genome-size " + s_genome_size_estimate + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
		
	
	#Смотрю, получился ли файл assembly.fasta. Его может не быть, если Flye не собрал ни одного дисджойнтига — в таком случае Flye прекращает работу преждевременно, не выдавая файла assembly.fasta. То, что Flye не выдал ни одного дисджойнтига, может быть связано с тем, что для ридов с большим количеством ошибок Mabs-flye попробовал очень маленький max_divergence. В случае, если файла assembly.fasta нет, я, даже не запуская скрипт calculate_AG.py, сразу считаю, что AG=0.
	if not os.path.isfile(s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta"):
		n_AG_for_point_2 = 0
	else:	
		
		#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
		os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + " --assembly " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta " + s_long_reads_option_for_calculate_AG + " " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe true --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

		#Беру AG, посчитанный скриптом calculate_AG.py
		if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt"):
			f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt", "r")
			s_line_1 = f_infile.readline()
			#AG is 487
			o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
			n_AG_for_point_2 = int(o_regular_expression_results.group(1))
		else:
			f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
			sys.exit()
	
	d_max_divergence_to_AG[n_max_divergence] = n_AG_for_point_2

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("AG for max_divergence " + str(n_max_divergence) + " is " + str(n_AG_for_point_2) + "\n\n")

	#Анализирую третью точку.
	n_number_of_the_point_under_analysis += 1
	n_max_divergence = n_point_3

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-flye started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. Max_divergence in this point is " + str(n_max_divergence) + "\n")

	#если пользователь не указывал размер генома
	if s_genome_size_estimate == "auto":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
	#если пользователь указал размер генома
	else:
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --genome-size " + s_genome_size_estimate + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)

	#Смотрю, получился ли файл assembly.fasta. Его может не быть, если Flye не собрал ни одного дисджойнтига — в таком случае Flye прекращает работу преждевременно, не выдавая файла assembly.fasta. То, что Flye не выдал ни одного дисджойнтига, может быть связано с тем, что для ридов с большим количеством ошибок Mabs-flye попробовал очень маленький max_divergence. В случае, если файла assembly.fasta нет, я, даже не запуская скрипт calculate_AG.py, сразу считаю, что AG=0.
	if not os.path.isfile(s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta"):
		n_AG_for_point_3 = 0
	else:
		#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
		os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + " --assembly " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta " + s_long_reads_option_for_calculate_AG + " " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe true --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

		#Беру AG, посчитанный скриптом calculate_AG.py
		if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt"):
			f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt", "r")
			s_line_1 = f_infile.readline()
			#AG is 487
			o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
			n_AG_for_point_3 = int(o_regular_expression_results.group(1))
		else:
			f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
			sys.exit()
		
	d_max_divergence_to_AG[n_max_divergence] = n_AG_for_point_3

	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("AG for max_divergence " + str(n_max_divergence) + " is " + str(n_AG_for_point_3) + "\n\n")

	#теперь последовательно выбираю остальные 8 точек методом золотого сечения и меряю AG для них.
	while n_number_of_the_point_under_analysis < 10: #"<", а не "<=", потому что увеличение номера точки здесь делается в начале цикла.
		n_number_of_the_point_under_analysis += 1
		
		#Смотрю, какая из двух центральных точек (вторая или третья) имеют меньшее значение AG. Если вторая имеет меньшее ли равное третьей, то выкидываю первую точку и сужаю интервал. Если третья имеет меньшее, чем вторая, то выкидываю четвёртую точку и сужаю интервал. При равных значениях выкидывыю левую, потому что равные значения могут быть из-за того, что в обеих точках AG = 0 из-за того, что при сборке по ридам с большим количеством ошибок обе центральные точки имели слишком маленький max_divergence, из-за чего Flye не смог собрать ни одного дисджойнтига и поэтому выдал пустой файл assembly.fasta.
		if n_AG_for_point_2 <= n_AG_for_point_3:
			n_point_1 = n_point_2
			n_point_2 = n_point_3
			#n_point_4 не меняется
			n_point_3 = round(10**(math.log10(n_point_4) - ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(math.log10(n_point_4) - math.log10(n_point_1))), 6)
			
			n_AG_for_point_1 = n_AG_for_point_2
			n_AG_for_point_2 = n_AG_for_point_3
			#n_AG_for_point_4 не меняется
			n_AG_for_point_3 = -100 #плейсхолдер. Всё равно это значение я сейчас посчитаю.
			
			#Анализирую третью точку.
			n_max_divergence = n_point_3

			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("Mabs-flye started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. Max_divergence in this point is " + str(n_max_divergence) + "\n")
			
			#если пользователь не указывал размер генома
			if s_genome_size_estimate == "auto":
				os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
			#если пользователь указал размер генома
			else:
				os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --genome-size " + s_genome_size_estimate + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
			
			#Смотрю, получился ли файл assembly.fasta. Его может не быть, если Flye не собрал ни одного дисджойнтига — в таком случае Flye прекращает работу преждевременно, не выдавая файла assembly.fasta. То, что Flye не выдал ни одного дисджойнтига, может быть связано с тем, что для ридов с большим количеством ошибок Mabs-flye попробовал очень маленький max_divergence. В случае, если файла assembly.fasta нет, я, даже не запуская скрипт calculate_AG.py, сразу считаю, что AG=0.
			if not os.path.isfile(s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta"):
				n_AG_for_point_3 = 0
			else:
				#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
				os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + " --assembly " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta " + s_long_reads_option_for_calculate_AG + " " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe true --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))

				#Беру AG, посчитанный скриптом calculate_AG.py
				if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt"):
					f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt", "r")
					s_line_1 = f_infile.readline()
					#AG is 487
					o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
					n_AG_for_point_3 = int(o_regular_expression_results.group(1))
				else:
					f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
					sys.exit()
			
			d_max_divergence_to_AG[n_max_divergence] = n_AG_for_point_3
			
			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("AG for max_divergence " + str(n_max_divergence) + " is " + str(n_AG_for_point_3) + "\n\n")
			
		elif n_AG_for_point_2 > n_AG_for_point_3:
			#n_point_1 не меняется
			n_point_4 = n_point_3
			n_point_3 = n_point_2
			n_point_2 = round(10**(math.log10(n_point_1) + ((math.sqrt(5) - 1) / (math.sqrt(5) + 1))*(math.log10(n_point_4) - math.log10(n_point_1))), 6)
			
			#n_AG_for_point_1 не меняется
			n_AG_for_point_4 = n_AG_for_point_3
			n_AG_for_point_3 = n_AG_for_point_2
			n_AG_for_point_2 = -100 #плейсхолдер. Всё равно это значение я сейчас посчитаю.
			
			#Анализирую вторую точку.
			n_max_divergence = n_point_2

			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("Mabs-flye started to analyze point " + str(n_number_of_the_point_under_analysis) + " of 10. Max_divergence in this point is " + str(n_max_divergence) + "\n")

			#если пользователь не указывал размер генома
			if s_genome_size_estimate == "auto":
				os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
			#если пользователь указал размер генома
			else:
				os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --out-dir " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + " --threads " + str(n_number_of_cpu_threads_to_use) + " --genome-size " + s_genome_size_estimate + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)			
			
			#Смотрю, получился ли файл assembly.fasta. Его может не быть, если Flye не собрал ни одного дисджойнтига — в таком случае Flye прекращает работу преждевременно, не выдавая файла assembly.fasta. То, что Flye не выдал ни одного дисджойнтига, может быть связано с тем, что для ридов с большим количеством ошибок Mabs-flye попробовал очень маленький max_divergence. В случае, если файла assembly.fasta нет, я, даже не запуская скрипт calculate_AG.py, сразу считаю, что AG=0.
			if not os.path.isfile(s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta"):
				n_AG_for_point_2 = 0
			else:
				#"--number_of_busco_orthogroups all" использую потому, что в папке BUSCO_dataset_to_use уже оставлены только те ортогруппы, которые нужно использовать.
				os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/calculate_AG.py --output_folder " + s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + " --assembly " + s_path_to_the_output_folder + "/Gene_assembly_for_max_divergence_" + str(n_max_divergence) + "/assembly.fasta " + s_long_reads_option_for_calculate_AG + " " + s_path_to_all_long_reads_that_correspond_to_busco_genes + " --number_of_busco_orthogroups all --local_busco_dataset " + s_path_to_the_output_folder + "/BUSCO_dataset_to_use --use_proovframe true --max_intron_length " + s_maximum_allowed_intron_length + " --threads " + str(n_number_of_cpu_threads_to_use))
				
				#Беру AG, посчитанный скриптом calculate_AG.py
				if os.path.isfile(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt"):
					f_infile = open(s_path_to_the_output_folder + "/AG_calculation_for_max_divergence_" + str(n_max_divergence) + "/AG.txt", "r")
					s_line_1 = f_infile.readline()
					#AG is 487
					o_regular_expression_results = re.search(r"AG is (\d+)", s_line_1)
					n_AG_for_point_2 = int(o_regular_expression_results.group(1))
				else:
					f_logs.write("Error. Couldn't calculate AG. See stderr and stdout for the reason why.")
			
			d_max_divergence_to_AG[n_max_divergence] = n_AG_for_point_2
			
			o_current_time_and_date = datetime.datetime.now()
			s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
			f_logs.write(s_current_time_and_date + "\n")
			f_logs.write("AG for max_divergence " + str(n_max_divergence) + " is " + str(n_AG_for_point_2) + "\n\n")
			
	
	#После того, как посчитал AG для всех 10 точек, я смотрю, какая из них имела наибольший AG. Тот max_divergence, который соответствует этой точке, я и считаю оптимальным для сборки Flye. Если две точки дают одинаковый AG, то, для определённости, выбираю ту из них, которая имеет больший max_divergence.
	
	n_max_divergence_that_provides_maximum_AG = -100
	n_maximum_AG = -100
	for n_max_divergence in d_max_divergence_to_AG:
		if d_max_divergence_to_AG[n_max_divergence] > n_maximum_AG:
			n_max_divergence_that_provides_maximum_AG = n_max_divergence
			n_maximum_AG = d_max_divergence_to_AG[n_max_divergence]
		
		if (d_max_divergence_to_AG[n_max_divergence] == n_maximum_AG) and (n_max_divergence > n_max_divergence_that_provides_maximum_AG):
			n_max_divergence_that_provides_maximum_AG = n_max_divergence
	
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("The optimal max_divergence is " + str(n_max_divergence_that_provides_maximum_AG) + ". When assembling only genes, it provides AG = " + str(n_maximum_AG) + ". Now Mabs-flye starts to assemble the genome using all reads and max_divergence = " + str(n_max_divergence_that_provides_maximum_AG) + "\n")
	
	#Теперь делаю сборку Flye по всем ридам, используя n_max_divergence_that_provides_maximum_AG.
	
	#если пользователь не указывал размер генома
	if s_genome_size_estimate == "auto":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_the_file_with_all_long_reads + " --out-dir " + s_path_to_the_output_folder + "/The_best_assembly --threads " + str(n_number_of_cpu_threads_to_use) + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence_that_provides_maximum_AG) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence_that_provides_maximum_AG) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
	#если пользователь указал размер генома
	else:
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Flye/bin/flye --nano-raw " + s_path_to_the_file_with_all_long_reads + " --out-dir " + s_path_to_the_output_folder + "/The_best_assembly --threads " + str(n_number_of_cpu_threads_to_use) + " --genome-size " + s_genome_size_estimate + " --no-alt-contigs --extra-params assemble_ovlp_divergence=" + str(n_max_divergence_that_provides_maximum_AG) + ",repeat_graph_ovlp_divergence=" + str(n_max_divergence_that_provides_maximum_AG) + ",assemble_divergence_relative=0 " + s_additional_flye_parameters)
	
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Mabs-flye finished. The contigs are in the file " + s_path_to_the_output_folder + "/The_best_assembly/assembly.fasta. Now I recommend to use a separate tool, for example HyPo (https://github.com/kensung-lab/hypo), to polish these contigs with accurate (HiFi, Illumina or MGI) reads.")





