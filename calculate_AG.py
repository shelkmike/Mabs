#!/usr/bin/env python3
# coding=utf-8

"""
calculate_AG, a part of the genome assembly suite Mabs. See https://github.com/shelkmike/Mabs

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
import statistics
import time
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

	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-map"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-map has not been found")
	
	if not os.path.isfile(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-fix"):
		l_unavailable_files_and_folders.append("The file " + s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-fix has not been found")


	#делаю парсинг аргументов командной строки. Можно было бы использовать argparse, но когда я делаю это без библиотек, то больше возможностей для того, чтобы сделать интерфейс таким, какой мне нравится.
		
	s_command_line = " ".join(sys.argv) #команда, которой запущен calculate_AG, в одну строку.
	s_command_line_reduced = s_command_line #то же, что s_command_line, но после того, как я распаршу какой-нибудь аргумент, я удалю его из этой строки. Если останется какой-то нераспарсенный аргумент, значит пользователь ввёл неизвестные calculate_AG аргументы, и нужно выдать ошибку.

	#инициализирую исходные значения переменных
	s_path_to_nanopore_reads = "" #путь к файлу с ридами Нанопора.
	s_path_to_pacbio_hifi_reads = "" #путь к файлу с ридами PacBio HiFi.
	s_path_to_pacbio_clr_reads = "" #путь к файлу с ридами PacBio CLR.
	s_busco_dataset_name_online = "" #название файла с базой данных BUSCO с сайта http://mikeshelk.site/Data/BUSCO_datasets/Latest/ (будет непустым только если пользователь использовал опцию "--download_busco_dataset") 
	s_path_to_a_local_busco_dataset = "" #путь к архивированному gzip файлу с датасетом BUSCO на диске или разархивированной папке с датасетом BUSCO на диске.
	n_number_of_cpu_threads_to_use = 10 #количество ядер, которые будет использовать calculate_AG.
	s_should_proovframe_be_used = "false" #нужно ли после шлифовки ридами использовать шлифовку Proovframe. Значение "true" или "false". 
	s_path_to_the_output_folder = "./AG_calculation_results" #путь к папке, куда calculate_AG пишет файлы.
	
	s_number_of_busco_orthogroups_to_use = "1000" #сколько ортогрупп BUSCO использовать. Это строка, содержащая или число, или слово "all", если нужно использовать все. Если пользователь укажет больше, чем есть в используемой базе данных BUSCO, то calculate_AG всё равно будет использовать все.
	s_maximum_allowed_intron_length = "from_BUSCO" #максимальная разрешённая длина интрона. По умолчанию, используется значение из файла dataset.cfg датасета BUSCO. Переменная начинается с "s_", потому что это строка. Ниже будет ещё переменная n_maximum_allowed_intron_length, которая число.

	s_version_of_calculate_AG = "2.12" #версия этой программы. Всегда равна версии Mabs. Поскольку эта программа нужна, в первую очередь, для Mabs, то когда я увеличиваю номер версии Mabs, то увеличивается и номер версии calculate_AG, и наоборот.

	l_errors_in_command_line = [] #список ошибок в командной строке. Если пользователь совершил много ошибок, то calculate_AG напишет про них все, а не только про первую встреченную.

	#если нет ни одного аргумента командной строки, или есть аргумент командной строки --help, то печатаю хелп
	if (len(sys.argv) == 1) or re.search(r"\s\-\-help", s_command_line):
		print("""calculate_AG, a program to assess how well genes were assembled in a genome assembly.

1) --assembly        Path to contigs or scaffolds to assess.
2) --nanopore_reads        Path to Nanopore reads.
3) --pacbio_clr_reads        Path to PacBio CLR reads, also known as "old PacBio" reads.
4) --pacbio_hifi_reads       Path to PacBio HiFi reads, also known as CCS reads.
5) --download_busco_dataset        Name of a file from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ . It should be the most taxonomically narrow dataset for your species. For example, for a human genome assembly, use "--download_busco_dataset primates_odb10.2021-02-19.tar.gz" and for a drosophila genome assembly use "--download_busco_dataset diptera_odb10.2020-08-05.tar.gz". Calculate_AG will download the respective file. This option is mutually exclusive with "--local_busco_dataset".
6) --threads        Number of CPU threads to be used by calculate_AG. The default value is 10.
7) --use_proovframe         Whether calculate_AG should polish sequences using Proovframe. "true" or "false". The default value is "false".
8) --output_folder        Output folder for calculate_AG results. The default is "AG_calculation_results".
9) --number_of_busco_orthogroups        How many BUSCO orthogroups should calculate_AG use. Should be either a positive integral value or "all" to use all orthogroups. The default value is 1000. 
10) --max_intron_length        Maximum allowed length of an intron. Should be either "from_BUSCO" to use a value from a BUSCO dataset, or a number, possibly ending with "k", "m" or "g". For example, 20k means 20 kilobases. The default is "from_BUSCO". Change --max_intron_length if you assemble a genome with unusually long introns.
11) --local_busco_dataset        Path to a local BUSCO dataset, manually pre-downloaded from http://mikeshelk.site/Data/BUSCO_datasets/Latest/ or http://busco-data.ezlab.org/v5/data/lineages/. Example: "--local_busco_dataset /home/test/Data/primates_odb10.2021-02-19.tar.gz". May be a .tar.gz or a decompressed folder. This option is mutually exclusive with "--download_busco_dataset".

Informational options:
12) --help        Print this help.
13) --version        Print the version of calculate_AG.

Example:
python3 calculate_AG.py --assembly contigs.fasta --nanopore_reads nanopore_reads.fastq --local_busco_dataset /mnt/lustre/username/Datasets/eudicots_odb10 --threads 40
		""")
		sys.exit()
		

	#смотрю, запросил ли пользователь версию calculate_AG
	if (len(sys.argv) == 1) or re.search(r"\s\-\-version", s_command_line):
		print("calculate_AG " + s_version_of_calculate_AG)
		sys.exit()

	#смотрю, дал ли пользователь последовательность сборки.
	o_regular_expression_results = re.search(r" --assembly (\S+)", s_command_line_reduced)
	if o_regular_expression_results:
		s_path_to_the_assembly = o_regular_expression_results.group(1)
		if not os.path.isfile(s_path_to_the_assembly):
			l_errors_in_command_line.append("The file with assembly " + s_path_to_the_assembly + " does not exist.")
		
		s_string_to_remove = re.escape(o_regular_expression_results.group(0))
		s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

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

	#проверяю, что пользователь дал хоть какие-то риды.
	if (s_path_to_nanopore_reads == "") and (s_path_to_pacbio_hifi_reads == "") and (s_path_to_pacbio_clr_reads == ""):
		l_errors_in_command_line.append("You need to provide at least one read set.")

	#создаю папку, в которой будут результаты calculate_AG
	o_regular_expression_results = re.search(r" --output_folder (\S+)", s_command_line_reduced)
	if o_regular_expression_results:
		s_path_to_the_output_folder = o_regular_expression_results.group(1)
		
		#Если в начале s_path_to_the_output_folder не стоит "." или "/", то, видимо, пользователь имеет в виду подпапку текущей папки, но не указал "./" в начале. В таком случае я добавлю "./" в начало, иначе могут быть проблемы. Не уверен насчёт mabs-flye и mabs-hifiasm, но у calculate_AG это вызывало проблемы.
		if (not re.search(r"^\.", s_path_to_the_output_folder)) and (not re.search(r"^\/", s_path_to_the_output_folder)):
			s_path_to_the_output_folder = "./" + s_path_to_the_output_folder
		
		s_string_to_remove = re.escape(o_regular_expression_results.group(0))
		s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
	
	#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка.  Не записываю эту ошибку в список l_errors_in_command_line , а сразу останавливаю работу, потому что если выходная папка уже существует, то в неё нельзя качать файлы BUSCO.
	if os.path.exists(s_path_to_the_output_folder):
		if len(os.listdir(s_path_to_the_output_folder)) != 0:
			print("calculate_AG has stopped because the output folder already exists and is not empty.")
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

	#смотрю, указал ли пользователь в командной строке, нужно ли использовать Proovframe.
	o_regular_expression_results = re.search(r" --use_proovframe (true|false|True|False|TRUE|FALSE)", s_command_line_reduced)
	if o_regular_expression_results:
		s_should_proovframe_be_used = o_regular_expression_results.group(1)
		#перевожу в нижний реестр.
		s_should_proovframe_be_used = s_should_proovframe_be_used.lower()
						
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

	#проверяю, не ввёл ли пользователь какие-то несуществующие опции. Это я определяю по тому, что после того, как я распарсил все команды, в строке s_command_line_reduced осталось что-то, кроме названия исполняемого файла calculate_AG и, возможно, пути перед ним.
	s_command_line_reduced = re.sub(r"^.*?calculate_AG(\.py)?\s*", "", s_command_line_reduced)
	if s_command_line_reduced != "":
		l_errors_in_command_line.append("You have provided some options which calculate_AG doesn't know: " + s_command_line_reduced)
		
	#проверяю, были ли недоступны какие-то программы, которые нужны Mabs, и были ли ошибки в командной строке. Если были какие-то из этих проблем, то пишу об этом и завершаю работу calculate_AG.
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
			print("There was an error in the command line of calculate_AG:")
			print(l_errors_in_command_line[0])
		#Если было больше одной ошибки.
		if len(l_errors_in_command_line) > 1:
			print("There were errors in the command line of calculate_AG:")
			n_error_number = 0 #порядковый номер ошибки. Считается от 1.
			for s_error_text in l_errors_in_command_line:
				n_error_number += 1
				print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
		
		#Печатаю пустую строку, как разделитель
		print("")
		
	if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
		#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме равно 1
		if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) == 1: 
			print("calculate_AG has stopped. Please, fix this error and restart calculate_AG.")
		#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме больше 1
		if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) > 1:
			print("calculate_AG has stopped. Please, fix these errors and restart calculate_AG.")
		
		sys.exit()

	################################
	#Со входными параметрами разобрался. Теперь, собственно, делаю работу.

	f_logs = open(s_path_to_the_output_folder + "/logs.txt","w",buffering=1)
	o_current_time_and_date = datetime.datetime.now()
	s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
	f_logs.write(s_current_time_and_date + "\n")
	f_logs.write("Started calculate_AG\n\n")

	f_logs.write("You have run calculate_AG of version " + s_version_of_calculate_AG + " with the following command: " + s_command_line + "\n\n")

	#сделаю специальный файл, в который в конце будет записана только строка вроде "AG is 1023".
	f_AG_calculation_results = open(s_path_to_the_output_folder + "/AG.txt", "w")
	
	#если пользователь сказал скачать файл с базой BUSCO или сам дал файл (но не папку), то разархивирую файл и меняю значение переменной s_path_to_a_local_busco_dataset с пути к файлу на путь к папке.
	if os.path.isfile(s_path_to_a_local_busco_dataset):
		s_busco_dataset_name = "" #в этой переменной название базы данных. Например, для файла eudicots_odb10.2020-09-10.tar.gz это будет "eudicots_odb10". После разархивирования файла, путь к которому содержится в переменной s_busco_dataset_archive, получается папка, название которой содержится в переменной s_busco_dataset_name.
		s_busco_dataset_name = re.sub(r"^.+\/", r"", s_path_to_a_local_busco_dataset) 
		s_busco_dataset_name = re.sub(r"\..+", r"", s_busco_dataset_name)
		os.system("tar -zxf " + s_path_to_a_local_busco_dataset + " --directory " + s_path_to_the_output_folder)
		
		s_path_to_a_local_busco_dataset = s_path_to_the_output_folder + "/" + s_busco_dataset_name
	
	#Оставляю из базы BUSCO только нужное количество (s_number_of_busco_orthogroups_to_use) ортогрупп — тех, которые имеют наиболее консервативные последовательности. Если пользователь указал использовать все ортогруппы, то calculate_AG использует все. Если пользователь указал больше ортогрупп, чем есть в этом наборе BUSCO, то calculate_AG использует все и пишет Warning в основной файл с логами.
	mabs_function_preprocess_busco_dataset.function_preprocess_busco_dataset(s_path_to_a_local_busco_dataset, s_number_of_busco_orthogroups_to_use, s_path_to_the_output_folder, f_logs)
	
	s_path_to_a_BUSCO_folder = s_path_to_the_output_folder + "/BUSCO_dataset_to_use/"

	#делаю ссылку на файл "ancestral", давая ему расширение .fasta. Затем делаю базу данных DIAMOND.
	#с помощью os.path.abspath() я получают абсолютный путь. Если он относительный, то это может создать проблемы в работоспособности мягкой ссылки.
	os.symlink(os.path.abspath(s_path_to_the_output_folder + "/BUSCO_dataset_to_use/ancestral"), s_path_to_the_output_folder + "/reference_busco_proteins.fasta")

	#Дальше я последовательности в сборке для простоты называю контигами, хотя они могут быть и скаффолдами.
		
	#если файл с контигами пустой, то сразу останавливаю выполнение calculate_AG, считая что AG=0. Иначе Metaeuk выдаст ошибку (если я правильно помню).
	n_size_of_the_file_with_contigs = os.stat(s_path_to_the_assembly).st_size
	if n_size_of_the_file_with_contigs == 0:
		f_logs.write("AG is 0")
		f_AG_calculation_results.write("AG is 0")
		sys.exit()

	#Делаю файл, который как входной, но в заголовках последовательностей все пробелы и табуляции заменены на "__". Иначе возникают проблемы, связанные с тем, что некоторые программы считают заголовком последовательности только то, что до первого пробельного символа.
	f_infile = open(s_path_to_the_assembly, "r")
	f_outfile = open(s_path_to_the_output_folder + "/assembly__without_whitespace_characters_in_titles.fasta", "w")
	for s_line in f_infile:
		if re.search(r"^>", s_line):
			s_title_line_without_whitespace_characters = re.sub(r"[ \t]", "__", s_line)
			f_outfile.write(s_title_line_without_whitespace_characters)
		else:
			f_outfile.write(s_line)
		
	f_infile.close()
	f_outfile.close()

	#Делаю шлифовку Proovframe, если пользователь попросил её сделать.
	if s_should_proovframe_be_used == "true":
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/DIAMOND/diamond makedb --in " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta -d " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta")
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-map --threads " + str(n_number_of_cpu_threads_to_use) + " --db " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta.dmnd --out " + s_path_to_the_output_folder + "/diamond_results_for_Proovframe.tsv " + s_path_to_the_output_folder + "/assembly__without_whitespace_characters_in_titles.fasta")
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Proovframe/bin/proovframe-fix --out " + s_path_to_the_output_folder + "/assembly_corrected_by_Proovframe.fasta " + s_path_to_the_output_folder + "/assembly__without_whitespace_characters_in_titles.fasta " + s_path_to_the_output_folder + "/diamond_results_for_Proovframe.tsv")
		
		s_path_to_the_assembly_after_polishing = s_path_to_the_output_folder + "/assembly_corrected_by_Proovframe.fasta"
	#Если пользователь не говорил сделать шлифовку Proovframe, то называю "сборкой после шлифовки" ту сборку, которую дал пользователь.
	else:
		s_path_to_the_assembly_after_polishing = s_path_to_the_output_folder + "/assembly__without_whitespace_characters_in_titles.fasta"

	#смотрю, какой в описании датасета BUSCO указан максимальный размер интрона. Это строка вида "max_intron=90000" в файле s_path_to_a_BUSCO_folder/dataset.cfg . Если пользователь сам указал эту величину через командную строку, то не беру значение из файла dataset.cfg.
	n_maximum_allowed_intron_length = -100 #плейсхолдер
	f_infile = open(s_path_to_a_BUSCO_folder + "/dataset.cfg", "r")
	for s_line in f_infile:
		if re.search(r"^max_intron\=(\d+)", s_line):
			o_regular_expression_results = re.search(r"^max_intron\=(\d+)", s_line)
			if s_maximum_allowed_intron_length == "from_BUSCO":
				n_maximum_allowed_intron_length = int(o_regular_expression_results.group(1))

	f_infile.close()


	#если пользователь сам указал конкретное значение максимального размера интрона
	if (s_maximum_allowed_intron_length != "from_BUSCO"):
		n_maximum_allowed_intron_length = int(s_maximum_allowed_intron_length)
	#если пользователь не указывал конкретное значение максимального размера интрона и n_maximum_allowed_intron_length всё ещё -100, значит это геном организма, у которого не бывает интронов, потому что в таком случае у BUSCO в файле .cfg нет информации про максимальный размер интрона. В таком случае я ставлю максимальный размер интронов 0.
	elif (n_maximum_allowed_intron_length == -100):
		n_maximum_allowed_intron_length = 0

		
	#делаю поиск генов BUSCO с помощью MetaEuk. MetaEuk пишет параллельно в количество файлов, равное количеству данных им потоков. Мои тесты на Макарыче показывают, что на Макарыче увеличение количества потоков выше 20 уже не приводят к увеличению скорости работы MetaEuk. А вот риск возникновения ошибки файловой системы увеличивается. Поэтому, на всякий случай, с MetaEuk я использую не более 20 потоков.
	n_number_of_cpu_threads_to_be_used_by_metaeuk = min(n_number_of_cpu_threads_to_use, 20)
	os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/MetaEuk/metaeuk easy-predict " + s_path_to_the_assembly_after_polishing + " " + s_path_to_the_output_folder + "/reference_busco_proteins.fasta " + s_path_to_the_output_folder + "/MetaEuk_results " + s_path_to_the_output_folder + "/MetaEuk_temporary_folder --threads " + str(n_number_of_cpu_threads_to_be_used_by_metaeuk) + " --max-intron " + str(n_maximum_allowed_intron_length) + " --overlap 1 --remove-tmp-files 1")

	#если MetaEuk вообще не выдал результатов, то считаю, что AG=0.
	if not os.path.exists(s_path_to_the_output_folder + "/MetaEuk_results.fas"):
		f_logs.write("AG is 0. Number of genes in single-copy orthogroups is 0. Number of genes in true multicopy orthogroups is 0. Number of genes in false multicopy orthogroups is 0.\n")
		f_AG_calculation_results.write("AG is 0")
		sys.exit()

	#если MetaEuk выдал результаты, но не нашёл ни одного гена, то считаю, что AG=0
	if os.path.getsize(s_path_to_the_output_folder + "/MetaEuk_results.fas") > 0:
		#теперь найденные Metaeuk гены анализирую с помощью HMMER, чтобы, соотнеся их с порогами BUSCO по bit score и по длине, понять, для каких ортогрупп сколько генов нашлось целиком, сколько фрагментировано, а для каких ортогрупп генов вообще не нашлось.
		s_path_to_MetaEuk_results = s_path_to_the_output_folder + "/MetaEuk_results.fas"
		s_path_to_the_file_with_BUSCO_scores_cutoff = s_path_to_a_BUSCO_folder + "/scores_cutoff"
		s_path_to_the_file_with_BUSCO_lengths_cutoff = s_path_to_a_BUSCO_folder + "/lengths_cutoff"

		#f_logs = open("logs.txt", "w", buffering = 1)

		#делаю словарь, в котором ключ это название ортогруппы, вроде 54443at71240, а значение это bit score cutoff, вроде 302.75.
		d_orthogroup_title_to_bit_score_cutoff = {}
		f_infile = open(s_path_to_the_file_with_BUSCO_scores_cutoff, "r")
		for s_line in f_infile:
			"""
			110468at71240   126.42
			84631at71240    378.98
			139803at71240   125.16
			"""
			
			o_regular_expression_results = re.search(r"^(\S+)\s+(\S+)", s_line)
			if o_regular_expression_results:
				s_orthogroup_title = o_regular_expression_results.group(1)
				n_bit_score_cutoff = float(o_regular_expression_results.group(2))
				d_orthogroup_title_to_bit_score_cutoff[s_orthogroup_title] = n_bit_score_cutoff
				
		f_infile.close()

		#делаю словари, в котором ключи это названия ортогрупп, вроде 54443at71240, а значение это, соответственно, средняя длина белка BUSCO и стандартное отклонение длин белков.
		d_orthogroup_title_to_the_average_BUSCO_protein_length = {}
		d_orthogroup_title_to_the_standard_deviation_of_BUSCO_protein_lengths = {}
		f_infile = open(s_path_to_the_file_with_BUSCO_lengths_cutoff, "r")
		for s_line in f_infile:
			"""
			Зачем нужен второй столбец, я не понимаю — насколько я вижу, в нём всегда ноль.
			110468at71240   0       25.5000316255   237
			84631at71240    0       31.4    314.0
			139803at71240   0       21.7    217.0
			"""
			o_regular_expression_results = re.search(r"^(\S+)\s+\S+\s+(\S+)\s+(\S+)", s_line)
			if o_regular_expression_results:
				s_orthogroup_title = o_regular_expression_results.group(1)
				n_length_standard_deviation = float(o_regular_expression_results.group(2))
				n_average_length = float(o_regular_expression_results.group(3))
				
				d_orthogroup_title_to_the_average_BUSCO_protein_length[s_orthogroup_title] = n_average_length
				d_orthogroup_title_to_the_standard_deviation_of_BUSCO_protein_lengths[s_orthogroup_title] = n_length_standard_deviation
				
		f_infile.close()
		
		#теперь померяю покрытие в контигах. При измерении покрытия использую только те риды, которые выравнивались к белкам BUSCO. Когда я меряю покрытие, я, для краткости, не выдаю позиции с нулевым покрытием. Это потому, что если геном достаточно большой, то, поскольку я выравниваю только те риды, которые выровнялись к белкам BUSCO, большинство позиций будут непокрытыми.
		#Если пользователь дал Mabs несколько наборов длинных ридов, то результаты их картирования сразу конкатенирую в один paf-файл mapping_to_calculate_coverage__before_simplification_of_split_reads.paf с помощью ">>".
		#У выходного файла окончание "__before_simplification_of_split_reads", потому что Minimap2 умеет делать split mapping, когда часть рида выравнивается в одно место, а часть в другое. При этом, split mapping может быть с частичным перекрытием фрагментов. Проблема в том, что я видел случаи, когда риды картировался как split mapping на ложно дуплицированные гены, которые были в разных контигах. При этом из-за существенного перекрытия картировавшихся частей у этих ридов было завышенное покрытие, из-за чего Mabs посчитал их истинно дуплицированными. Поэтому ниже я буду изменять файл mapping_to_calculate_coverage__before_simplification_of_split_reads.paf , превращая его в файл mapping_to_calculate_coverage.paf , при этом для каждого рида будет оставляться только одно картирование.
		
		#если файл mapping_to_calculate_coverage__before_simplification_of_split_reads.paf существует, то удаляю его. А то, если я случайно при запуске calculate_AG не удалил прошлую папку с результатами, то будут проблемы из-за дозаписи в конец существующего файла. Впрочем, думаю что при неудалении старой папки возможны и какие-то другие проблемы. Потом, возможно, нужно будет сделать у calculate_AG какое-нибудь предохранение на случай, если пользователь даёт в качестве выходной непустую папку.
		if os.path.isfile(s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf"):
			os.remove(s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf")
		
		#если пользователь дал риды Нанопора.
		if (s_path_to_nanopore_reads != ""):
			os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2 -x map-ont -I 1000G -t " + str(n_number_of_cpu_threads_to_use) + " " + s_path_to_the_assembly_after_polishing + " " + s_path_to_nanopore_reads + " >>" + s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf")

		#если пользователь дал риды CLR.
		if (s_path_to_pacbio_clr_reads != ""):
			os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2 -x map-pb -I 1000G -t " + str(n_number_of_cpu_threads_to_use) + " " + s_path_to_the_assembly_after_polishing + " " + s_path_to_pacbio_clr_reads + " >>" + s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf")
		
		#если пользователь дал риды HiFi.
		if (s_path_to_pacbio_hifi_reads != ""):
			os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Minimap2/minimap2 -x map-hifi -I 1000G -t " + str(n_number_of_cpu_threads_to_use) + " " + s_path_to_the_assembly_after_polishing + " " + s_path_to_pacbio_hifi_reads + " >>" + s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf")	
				
		#Теперь нужно изменить paf-файл так, чтобы для каждого рида осталось только одно картирование. См. пояснение выше, после фразы "У выходного файла окончание "__before_simplification_of_split_reads", потому что...".
		f_infile = open(s_path_to_the_output_folder + "/mapping_to_calculate_coverage__before_simplification_of_split_reads.paf", "r")
		f_outfile = open(s_path_to_the_output_folder + "/mapping_to_calculate_coverage.paf", "w")

		s_title_of_the_previous_read = "" #заголовок предыдущего рида.

		for s_line in f_infile:
			#m140928_184123_42139_c100719602550000001823155305141590_s1_p0/66/0_11284        11284   291     11185   -       utg000045l      2349263 2300607 2312221 2672    11710   60      tp:A:P  cm:i:164        s1:i:2374       s2:i:0  dv:f:0.0749     rl:i:52
			#m140928_184123_42139_c100719602550000001823155305141590_s1_p0/93/0_17242        17242   29      17200   +       utg000004l      1653227 994564  1010965 5961    17430   60      tp:A:P  cm:i:337        s1:i:5568       s2:i:135        dv:f:0.0522     rl:i:623

			l_line_split = s_line.split("\t")
			if re.search(r"^\d+$", l_line_split[1]): #если я вижу во 2-м столбце число, значит это строка, описывающая картирование рида на контиг.
				n_read_title = l_line_split[0]
				
				if n_read_title != s_title_of_the_previous_read: #если в этой строке описывается картирование не того же рида, картирование которого описывалось в прошлой строке.
					f_outfile.write(s_line)	
				s_title_of_the_previous_read = n_read_title
				
		f_infile.close()
		f_outfile.close()
		
		#конвертирую paf в bed. согласно описанию формата paf в восьмом столбце "Target start on original strand (0-based)", а в девятом столбце "Target end on original strand (0-based)".
		"""
		ERR5530736.228149       16542   176     16529   +       utg000209l      489714  242488  259229  7174    16836   60      tp:A:P  cm:i:593        s1:i:7085       s2:i:97 dv:f:0.0178     rl:i:477
		ERR5530736.120494       19723   237     19619   +       utg000165l      261240  133700  153928  2716    20331   60      tp:A:P  cm:i:152        s1:i:2480       s2:i:165        dv:f:0.0511     rl:i:296
		ERR5530736.132409       3421    239     3318    +       utg000063l      324047  218178  221319  445     3152    60      tp:A:P  cm:i:28 s1:i:418        s2:i:89 dv:f:0.0453     rl:i:0
		ERR5530736.682110       7703    174     7633    -       utg000275l      77109   59845   67405   2838    7633    60      tp:A:P  cm:i:233        s1:i:2764       s2:i:1223       dv:f:0.0176     rl:i:102
		"""

		f_infile = open(s_path_to_the_output_folder + "/mapping_to_calculate_coverage.paf", "r")
		f_outfile = open(s_path_to_the_output_folder + "/mapping_to_calculate_coverage.bed", "w")

		for s_line in f_infile:
			if not re.search(r"^\s*$", s_line):
				l_line_split = s_line.split("\t")
				s_contig_title = l_line_split[5]
				#в формате paf первая из двух координат всегда меньше второй
				n_match_leftmost_coordinate = int(l_line_split[7])
				n_match_rightmost_coordinate = int(l_line_split[8])
				f_outfile.write(s_contig_title + "\t" + str(n_match_leftmost_coordinate) + "\t" + str(n_match_rightmost_coordinate + 1) + "\n")

		f_infile.close()
		f_outfile.close()

		#сортирую bed-файл с ридами. Сортировка рекомедуется для ускорения поиска перекрытий на https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html .
		os.system("sort -k 1,1 -k2,2n " + s_path_to_the_output_folder + "/mapping_to_calculate_coverage.bed > " + s_path_to_the_output_folder + "/mapping_to_calculate_coverage.sorted.bed")
		
		#Для экономии времени, я при расчёте покрытия буду использовать не весь bed-файл с картированиями ридов, а только его пересечение с bed-файлом экзонов, найденных Metaeuk. Потому что покрытие не в экзонах меня всё равно не интересует.
		#Делаю bed-файл с координатами экзонов, найденных Metaeuk. Тут не только те экзоны, которые принадлежат генам, прошедшим фильтрацию по длине и bit score BUSCO, а все экзоны всех генов, выданных Metaeuk. Это для простоты. Возможно, для скорости потом стоит сделать, чтобы использовались только те экзоны, которые принадлежат к Complete генам, согласно критериям по длине и bit score.
		f_infile = open(s_path_to_the_output_folder + "/MetaEuk_results.fas", "r")
		f_outfile = open(s_path_to_the_output_folder + "/all_exons_found_by_MetaEuk.bed", "w")
		
		for s_line in f_infile:
			"""
			#Заголовки имеют вид 192750at33090|utg002705l|-|115|1.18e-28|3|19694|20872|20872[20872]:20765[20765]:108[108]|20515[20515]:20336[20336]:180[180]|19924[19924]:19694[19694]:231[231]
			#На https://github.com/soedinglab/metaeuk написано, что 19694|20872 это координаты начала и конца гена. 20872[20872]:20765[20765] это координаты начала и конца экзона. Если два экзона перекрываются (metaeuk почему-то позволяет такое), то значения за квадратными скобками будут включать и перекрывающиеся нуклеотиды, а внутри квадратных скобок - нет. Поэтому я буду брать значения внутри квадратных скобок.
			#Metaeuk записывает координаты экзонов, считая их от 0, а не от 1: на https://github.com/soedinglab/metaeuk написано "coord refers to the coordination on the contig (first base has coordinate 0)".
			"""
			o_regular_expression_results = re.search(r"^>(.*?)\|([^\|]+)\|(\+|\-)\|.*?\|.*?\|.*?\|(\d+)\|(\d+)\|(.+)", s_line)

			if o_regular_expression_results:
				s_orthogroup_title = o_regular_expression_results.group(1)
				s_contig_title = o_regular_expression_results.group(2)
				s_chain = o_regular_expression_results.group(3) #цепь, на которой лежит ген. "+" или "-".
				n_leftmost_coordinate_of_the_gene = int(o_regular_expression_results.group(4)) + 1 #прибавляю единицу, потому что Metaeuk выдаёт координаты zero-based.
				n_rightmost_coordinate_of_the_gene = int(o_regular_expression_results.group(5))
				s_string_with_exons = o_regular_expression_results.group(6)
				
				#иду по всем координатам экзонов. Координаты экзона содержатся в подстроке вида 20872[20872]:20765[20765]:108[108]
				while re.search(r"\d+\[(\d+)\]\:\d+\[(\d+)\]\:\d+\[\d+\]", s_string_with_exons):
					o_regular_expression_results_2 = re.search(r"\d+\[(\d+)\]\:\d+\[(\d+)\]\:\d+\[\d+\]", s_string_with_exons)
					n_first_exon_coordinate = int(o_regular_expression_results_2.group(1))
					n_second_exon_coordinate = int(o_regular_expression_results_2.group(2))
					
					#если ген обратно-комплементарный, то первой координатой была записана бОльшая. Делаю первой координатой меньшую.
					if n_first_exon_coordinate > n_second_exon_coordinate:
						n_temp = n_second_exon_coordinate
						n_second_exon_coordinate = n_first_exon_coordinate
						n_first_exon_coordinate = n_temp
					
					f_outfile.write(s_contig_title + "\t" + str(n_first_exon_coordinate) + "\t" + str(n_second_exon_coordinate) + "\n")
					
					n_first_exon_coordinate += 1 #прибавляю единицу, потому что Metaeuk выдаёт координаты zero-based.
										
					#удаляю упоминание об этом экзоне из строки, чтобы можно было начать рассматривать новый.
					#f_logs.write("analyzed exon " + o_regular_expression_results_2.group(0) + "\n")
					s_exon_information_with_masked_metacharacters = re.escape(o_regular_expression_results_2.group(0)) #s_exon_information_with_masked_metacharacters это как o_regular_expression_results_2.group(0) , но все метасимволы замаскированы. Нужно, чтобы правильно прошло удаление этой подстроки из s_string_with_exons с помощью re.sub
					s_string_with_exons = re.sub(s_exon_information_with_masked_metacharacters, "", s_string_with_exons)
		
		f_infile.close()
		f_outfile.close()
		
		#сортирую bed-файл с координатами экзонов. Сортировка рекомедуется для ускорения поиска перекрытий на https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html .
		os.system("sort -k 1,1 -k2,2n " + s_path_to_the_output_folder + "/all_exons_found_by_MetaEuk.bed > " + s_path_to_the_output_folder + "/all_exons_found_by_MetaEuk.sorted.bed")
		
		#теперь ищу перекрытия между экзонами и ридами.
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Bedtools/bedtools intersect -sorted -a " + s_path_to_the_output_folder + "/mapping_to_calculate_coverage.sorted.bed -b " + s_path_to_the_output_folder + "/all_exons_found_by_MetaEuk.sorted.bed >" + s_path_to_the_output_folder + "/parts_of_reads_that_cover_exons.bed")
		
		#сортирую как написано на https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
		os.system("sort -k 1,1 " + s_path_to_the_output_folder + "/parts_of_reads_that_cover_exons.bed >" + s_path_to_the_output_folder + "/parts_of_reads_that_cover_exons.sorted.bed")

		#генерирую специальный геномный файл, который нужен Bedtools. В нём для каждого контига идёт заголовок, а потом, через табуляцию, длина этого контига.
		f_outfile = open(s_path_to_the_output_folder + "/regions.sizes", "w")
		
		#сначала сделаю словарь, где ключ это заголовок контига без ">", а значение — длина контига.
		d_contig_title_to_contig_length = {} 
		f_infile = open(s_path_to_the_assembly_after_polishing, "r")
		
		s_contig_title = ""
		for s_line in f_infile:
			o_regular_expression_results = re.search(r"^>(.+)", s_line)
			if o_regular_expression_results:
				s_contig_title = o_regular_expression_results.group(1)
			elif re.search(r"^(.+)", s_line): #если это не строка с заголовком, то считаю, что строка с последовательностью
				o_regular_expression_results = re.search(r"^(.+)",s_line)
				
				#если для этого контига значение длины ещё не инициализировано, то инициализирую его
				if s_contig_title not in d_contig_title_to_contig_length:
					d_contig_title_to_contig_length[s_contig_title] = 0

				s_sequence_in_this_string = o_regular_expression_results.group(1)
				#удаляю пробельные символы
				s_sequence_in_this_string = re.sub(r"\s", "", s_sequence_in_this_string)
				d_contig_title_to_contig_length[s_contig_title] += len(s_sequence_in_this_string)
		f_infile.close()
		
		for s_contig_title in d_contig_title_to_contig_length:
			f_outfile.write(s_contig_title + "\t" + str(d_contig_title_to_contig_length[s_contig_title]) + "\n")
		
		f_outfile.close()
		
		#считаю распределение покрытия
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/Bedtools/bedtools genomecov -dz -i " + s_path_to_the_output_folder + "/parts_of_reads_that_cover_exons.sorted.bed -g " + s_path_to_the_output_folder + "/regions.sizes >" + s_path_to_the_output_folder + "/list_of_coverages.txt")
		
		#теперь записываю покрытие всех позиций, покрытых хотя бы одним ридом, в двойной словарь, где первый ключ это заголовок контига, а второй это позиция, например ["utg000098l"]["763"] , а значение это покрытие длинными ридами из файла long_reads_that_have_matches_to_BUSCO_proteins.fastq в этой позиции. Позиции отсчитываются от 1.
		dd_contig_title_and_position_to_coverage = {}
		f_infile = open(s_path_to_the_output_folder + "/list_of_coverages.txt", "r")
		for s_line in f_infile:
			#Координаты в этом списке 0-based.
			"""
			utg000107l_163374-203511        17321   51
			utg000107l_163374-203511        17322   51
			utg000107l_163374-203511        17323   50
			utg000107l_163374-203511        17324   51
			utg000107l_163374-203511        17325   53
			"""
			if re.search(r"^(\S+)\s+(\S+)\s+(\S+)", s_line):
				o_regular_expression_results = re.search(r"^(\S+)\s+(\S+)\s+(\S+)", s_line)
				s_region_title = o_regular_expression_results.group(1)
				n_position = int(o_regular_expression_results.group(2)) + 1
				n_coverage = int(o_regular_expression_results.group(3))
				
				#инициализирую словарь для этого контига, если он ещё не инициализирован
				if s_region_title not in dd_contig_title_and_position_to_coverage:
					dd_contig_title_and_position_to_coverage[s_region_title] = {}
				
				dd_contig_title_and_position_to_coverage[s_region_title][n_position] = n_coverage
		
		f_infile.close()
		
		#os.system("rm " + s_path_to_the_output_folder + "/list_of_coverages.txt") #чтобы освободить место
		
		
		#теперь я делаю выравнивание всех марковских моделей к белкам, найденным Metaeuk. По результатам hmmsearch я определяю, сколько для каждой ортогруппы complete генов, сколько fragmented и сколько missing.
		
		#Параллельно, составлю такие два словаря списков:
		dl_orthogroup_title_to_the_list_of_its_genes = {} #ключ - название ортогруппы, вроде "192750at33090". Значение — список с описаниями генов. Описание включает название контига, которому принадлежит ген, самую левую координату гена и самую правую координату гена, вроде "utg002705l:19694-20872". В этом словаре списков только гены, которые собрались целиком, то есть не фрагментированно.
		#ВАЖНО! В строках вида "utg002705l:19694-20872" я использую one-based координаты, а не zero-based. Соответственно, поскольку metaeuk, paf и bed используют zero-based формат, то мне нужно периодически конвертировать из одного в другое.
		
		#это словарь списков, в котором ключ это название контига и координаты гена внутри конига, вроде "utg002705l:19694-20872", а значение это список покрытий всех позиций всех экзонов этого гена.
		dl_gene_description_to_the_list_of_coverages_in_its_exons = {}
		
				
		#Делаю поиск hmmsearch с параметрами по умолчанию — так же, как делает BUSCO. По умолчанию порог по e-value 10. "-o /dev/null" использую, чтобы hmmsearch не писал кучу лишней информации в стандартный вывод.
		os.system(s_path_to_the_folder_where_Mabs_lies + "/Additional/HMMER/src/hmmsearch --domtblout " +  s_path_to_the_output_folder + "/hmmsearch_results.txt -o /dev/null --cpu " + str(n_number_of_cpu_threads_to_use) + " " + s_path_to_a_BUSCO_folder + "/concatenated_profile_HMMs_of_orthogroups.hmm " + s_path_to_MetaEuk_results)
			
		f_infile = open(s_path_to_the_output_folder + "/hmmsearch_results.txt", "r")
		
		dl_orthogroup_title_to_the_list_of_targets_I_have_already_seen_in_this_file = {} #словарь списков, в котором ключ это название ортогруппы, вроде 164894at71240, а список содержит те мишени, которые я уже видел в этом файле. Дело в том, что hmmsearch делает локальное выравнивание и для выравнивания одного марковского профиля к одной таргетной последовательности может быть несколько строк. Впрочем, в каждой строке есть информация и по длине всей таргетной последовательности и по тому bit score, который получается по всем локальным выравниваниям этой таргетной последовательности вместе взятой, поэтому без разницы, на какую строку для данной ортогруппы смотреть.
		
		for s_line in f_infile:
			"""
			#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
			# target name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
			#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
			164894at71240|utg000016l|+|359|2.123e-102|5|2630724|2631813|2630724[2630724]:2631002[2631002]:279[279]|2631092[2631098]:2631160[2631160]:69[63]|2631237[2631246]:2631305[2631305]:69[60]|2631498[2631498]:2631647[2631647]:150[150]|2631736[2631739]:2631813[2631813]:78[75]                                                                                                                                                                                                                                                                                                                                                   -            209 164894at71240        -            165   4.7e-73  241.8   0.1   1   2     0.056        55    0.9   0.0     1    24     1    24     1    45 0.80 -
			164894at71240|utg000016l|+|359|2.123e-102|5|2630724|2631813|2630724[2630724]:2631002[2631002]:279[279]|2631092[2631098]:2631160[2631160]:69[63]|2631237[2631246]:2631305[2631305]:69[60]|2631498[2631498]:2631647[2631647]:150[150]|2631736[2631739]:2631813[2631813]:78[75]                                                                                                                                                                                                                                                                                                                                                   -            209 164894at71240        -            165   4.7e-73  241.8   0.1   2   2   3.3e-75   3.2e-72  239.1   0.0    13   160    52   205    41   209 0.88 -
			47967at71240|utg000009l|-|1175|0|12|11472827|11476317|11476317[11476317]:11475961[11475961]:357[357]|11475875[11475875]:11475792[11475792]:84[84]|11475703[11475703]:11475620[11475620]:84[84]|11475528[11475522]:11475469[11475469]:60[54]|11475378[11475378]:11475325[11475325]:54[54]|11475224[11475218]:11475078[11475078]:147[141]|11474968[11474965]:11474867[11474867]:102[99]|11474641[11474638]:11474525[11474525]:117[114]|11474425[11474422]:11474294[11474294]:132[129]|11474118[11474118]:11473867[11473867]:252[252]|11473752[11473752]:11472958[11472958]:795[795]|11472874[11472868]:11472827[11472827]:48[42] -            735 164894at71240        -            165   5.4e-06   23.7   0.4   1   2   1.8e-08   1.7e-05   22.0   0.0    82   147   134   201   128   204 0.87 -
			47967at71240|utg000009l|-|1175|0|12|11472827|11476317|11476317[11476317]:11475961[11475961]:357[357]|11475875[11475875]:11475792[11475792]:84[84]|11475703[11475703]:11475620[11475620]:84[84]|11475528[11475522]:11475469[11475469]:60[54]|11475378[11475378]:11475325[11475325]:54[54]|11475224[11475218]:11475078[11475078]:147[141]|11474968[11474965]:11474867[11474867]:102[99]|11474641[11474638]:11474525[11474525]:117[114]|11474425[11474422]:11474294[11474294]:132[129]|11474118[11474118]:11473867[11473867]:252[252]|11473752[11473752]:11472958[11472958]:795[795]|11472874[11472868]:11472827[11472827]:48[42] -            735 164894at71240        -            165   5.4e-06   23.7   0.4   2   2      0.43   4.2e+02   -2.0   0.1    93   139   575   622   572   644 0.66 -
			100085at71240|utg000009l|-|405|3.016e-116|1|7632185|7632937|7632937[7632937]:7632185[7632185]:753[753]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         -            251 164894at71240        -            165     0.014   12.6   0.0   1   1   2.1e-05      0.02   12.1   0.0   104   139     2    37     1    49 0.89 -
			141127at71240|utg000021l|-|502|1.904e-145|9|1743371|1746261|1746261[1746261]:1745905[1745905]:357[357]|1745323[1745323]:1745258[1745258]:66[66]|1745146[1745146]:1745081[1745081]:66[66]|1744743[1744728]:1744645[1744645]:99[84]|1744128[1744122]:1744006[1744006]:123[117]|1743894[1743891]:1743844[1743844]:51[48]|1743755[1743752]:1743702[1743702]:54[51]|1743619[1743613]:1743551[1743551]:69[63]|1743406[1743406]:1743371[1743371]:36[36]                                                                                                                                                                               -            296 164894at71240        -            165      0.04   11.1   0.0   1   1   0.00013      0.12    9.5   0.0   104   140   247   283   190   287 0.75 -
			#
			# Program:         hmmsearch
			# Version:         3.1b2 (February 2015)
			# Pipeline mode:   SEARCH
			# Query file:      /mnt/lustre/shelkmike/Work/Data/Protein_and_nucleotide_datasets/BUSCO/V5/eudicots_odb10/hmms/164894at71240.hmm
			# Target file:     ./MetaEuk_results.fas
			# Option settings: /mnt/lustre/shelkmike/Work/Tools/HMMER/hmmer-3.1b2/Installed_here_by_mikeshelk/bin/hmmsearch --domtblout ./Hmmsearch_results/hmmsearch_results_for164894at71240.txt --cpu 22 /mnt/lustre/shelkmike/Work/Data/Protein_and_nucleotide_datasets/BUSCO/V5/eudicots_odb10/hmms/164894at71240.hmm ./MetaEuk_results.fas 
			# Current dir:     /mnt/lustre/shelkmike/Work/Mine/My_tools/Mabs/1.8/Test_alternative_variant_35/Test_MetaEuk
			# Date:            Fri Feb 11 14:25:03 2022
			# [ok]
			"""
			
			#Заголовки имеют вид 192750at33090|utg002705l|-|115|1.18e-28|3|19694|20872|20872[20872]:20765[20765]:108[108]|20515[20515]:20336[20336]:180[180]|19924[19924]:19694[19694]:231[231]
			#На https://github.com/soedinglab/metaeuk написано, что 19694|20872 это координаты начала и конца гена. 20872[20872]:20765[20765] это координаты начала и конца экзона. Если два экзона перекрываются (metaeuk почему-то позволяет такое), то значения за квадратными скобками будут включать и перекрывающиеся нуклеотиды, а внутри квадратных скобок - нет. Поэтому я буду брать значения внутри квадратных скобок.
			#Metaeuk записывает координаты экзонов, считая их от 0, а не от 1: на https://github.com/soedinglab/metaeuk написано "coord refers to the coordination on the contig (first base has coordinate 0)".
			
			#если это не пустая строка и не строка с комментарием
			if (not re.search("^#", s_line)) and (re.search("\S", s_line)):
				s_line = re.sub(r"[\r\n]+$", "", s_line)
				l_line_split = s_line.split() #когда split используется без аргументов, Питон разбивает строку по группам идущих подряд пробельных символов.
				s_target_name = l_line_split[0]
				n_target_length = int(l_line_split[2])
				s_orthogroup_title = l_line_split[3]
				n_bit_score = float(l_line_split[7])
				
				if s_orthogroup_title not in dl_orthogroup_title_to_the_list_of_targets_I_have_already_seen_in_this_file:
					dl_orthogroup_title_to_the_list_of_targets_I_have_already_seen_in_this_file[s_orthogroup_title] = []
				
				if s_target_name not in dl_orthogroup_title_to_the_list_of_targets_I_have_already_seen_in_this_file[s_orthogroup_title]:
					dl_orthogroup_title_to_the_list_of_targets_I_have_already_seen_in_this_file[s_orthogroup_title].append(s_target_name)
					
					if n_bit_score >= d_orthogroup_title_to_bit_score_cutoff[s_orthogroup_title]:
						n_z_value = (n_target_length - d_orthogroup_title_to_the_average_BUSCO_protein_length[s_orthogroup_title]) / d_orthogroup_title_to_the_standard_deviation_of_BUSCO_protein_lengths[s_orthogroup_title]
						
						if n_z_value >= -2:													
							o_regular_expression_results = re.search(r"^(.*?)\|([^\|]+)\|(\+|\-)\|.*?\|.*?\|.*?\|(\d+)\|(\d+)\|(.+)", s_target_name)

							if o_regular_expression_results:
								s_orthogroup_title = o_regular_expression_results.group(1)
								s_contig_title = o_regular_expression_results.group(2)
								s_chain = o_regular_expression_results.group(3) #цепь, на которой лежит ген. "+" или "-".
								n_leftmost_coordinate_of_the_gene = int(o_regular_expression_results.group(4)) + 1 #прибавляю единицу, потому что Metaeuk выдаёт координаты zero-based.
								n_rightmost_coordinate_of_the_gene = int(o_regular_expression_results.group(5))
								s_string_with_exons = o_regular_expression_results.group(6)
								
								s_gene_description = s_contig_title + ":" + str(n_leftmost_coordinate_of_the_gene) + "-" + str(n_rightmost_coordinate_of_the_gene)
								
								if s_orthogroup_title not in dl_orthogroup_title_to_the_list_of_its_genes:
									dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title] = []
								dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title].append(s_gene_description)
								
								#f_logs.write("Started to analyze the coverage in exons of " + s_gene_description + "\n")
								
								if s_gene_description not in dl_gene_description_to_the_list_of_coverages_in_its_exons:
									dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description] = []
								
								#иду по всем координатам экзонов. Координаты экзона содержатся в подстроке вида 20872[20872]:20765[20765]:108[108]
								while re.search(r"\d+\[(\d+)\]\:\d+\[(\d+)\]\:\d+\[\d+\]", s_string_with_exons):
									o_regular_expression_results_2 = re.search(r"\d+\[(\d+)\]\:\d+\[(\d+)\]\:\d+\[\d+\]", s_string_with_exons)
									n_first_exon_coordinate = int(o_regular_expression_results_2.group(1))
									n_second_exon_coordinate = int(o_regular_expression_results_2.group(2))
									
									#если ген обратно-комплементарный, то первой координатой была записана бОльшая. Делаю первой координатой меньшую.
									if n_first_exon_coordinate > n_second_exon_coordinate:
										n_temp = n_second_exon_coordinate
										n_second_exon_coordinate = n_first_exon_coordinate
										n_first_exon_coordinate = n_temp
									
									n_first_exon_coordinate += 1 #прибавляю единицу, потому что Metaeuk выдаёт координаты zero-based.
									
									l_coverages_in_this_exon = [] #список покрытий в этом экзоне
									
									for n_position in range(n_first_exon_coordinate, n_second_exon_coordinate + 1):
										
										#проверяю, нет ли такого, что в этом контиге не было вообще ни одной покрытой позиции. Тогда dd_contig_title_and_position_to_coverage будет неинициализирован для s_contig_title. В таком случае я считаю покрытие в этой позиции равным нулю.
										if s_contig_title not in dd_contig_title_and_position_to_coverage:
											n_coverage = 0
										else:
											#поскольку в двойной словарь dd_contig_title_and_position_to_coverage я записывал покрытие только тех позиций, покрытие которых было ненулевым, то сейчас нужно проверить, есть ли эта позиция в этом двойном словаре.
											if n_position in dd_contig_title_and_position_to_coverage[s_contig_title]:
												n_coverage = dd_contig_title_and_position_to_coverage[s_contig_title][n_position]
											else:
												n_coverage = 0
											
										l_coverages_in_this_exon.append(n_coverage)
										
									#добавляю список покрытий этого экзона к словарю списков, который содержит списки покрытий для каждого гена
									dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description] += l_coverages_in_this_exon									
									#удаляю упоминание об этом экзоне из строки, чтобы можно было начать рассматривать новый.
									#f_logs.write("analyzed exon " + o_regular_expression_results_2.group(0) + "\n")
									s_exon_information_with_masked_metacharacters = re.escape(o_regular_expression_results_2.group(0)) #s_exon_information_with_masked_metacharacters это как o_regular_expression_results_2.group(0) , но все метасимволы замаскированы. Нужно, чтобы правильно прошло удаление этой подстроки из s_string_with_exons с помощью re.sub
									s_string_with_exons = re.sub(s_exon_information_with_masked_metacharacters, "", s_string_with_exons)
						
						#если ген присутствует, но фрагментирован, то я информацию о нём никак не использую.
						else:
							pass

		f_infile.close()
		
		#если ни один из найденных генов не является complete, то я считаю, что AG равен 0, и перехожу к вычислению N50. То, что ни один из генов не complete, я определяю по тому, что в dl_orthogroup_title_to_the_list_of_its_genes нет информации ни для одной ортогруппы.
		if len(dl_orthogroup_title_to_the_list_of_its_genes) == 0:
			n_AG = 0
			n_number_of_single_copy_genes_found_in_the_assembly = 0
			n_number_of_true_multicopy_genes = 0
			n_number_of_false_multicopy_genes = 0
		else:
			#Теперь посчитаю медианное покрытие по всем экзонам всех однокопийных ортогрупп. Для этого сначала составлю список, в котором будут конкатенированы все списки покрытий по однокопийным ортогруппам. Параллельно, посчитаю количество однокопийных ортогрупп.
			#Кроме этого, параллельно я сделаю словарь, в котором в соответствие каждому гену (строки вроде "utg002705l:19694-20872") послевлено медианное покрытие его экзонов.
			d_gene_description_to_the_median_coverage_in_its_exons = {}
			l_coverages_in_exons_of_single_copy_genes = []
			l_coverages_in_exons_of_genes = [] #как l_coverages_in_exons_of_single_copy_genes, но считается не по однокопийным генам, а по всем генам.
			n_number_of_single_copy_genes_found_in_the_assembly = 0
			for s_orthogroup_title in dl_orthogroup_title_to_the_list_of_its_genes:
				
				#если для ортогруппы найден всего один ген.
				if len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title]) == 1:
					s_gene_description = dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title][0]
					
					d_gene_description_to_the_median_coverage_in_its_exons[s_gene_description] = statistics.median(dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description])
					
					l_coverages_in_exons_of_single_copy_genes += dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description]
					n_number_of_single_copy_genes_found_in_the_assembly += 1
					
					l_coverages_in_exons_of_genes += dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description]
					
					f_logs.write("For the gene " + s_gene_description + " from a single-copy orthogroup " + s_orthogroup_title + ", the median coverage is " + str(d_gene_description_to_the_median_coverage_in_its_exons[s_gene_description]) + ". It was calculated using " + str(len(dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description])) + " positions.\n")
				
				#если для ортогруппы найдено больше одного гена.
				if len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title]) > 1:
					for s_gene_description in dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title]:
						d_gene_description_to_the_median_coverage_in_its_exons[s_gene_description] = statistics.median(dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description])
						
						l_coverages_in_exons_of_genes += dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description]
						
						f_logs.write("For the gene " + s_gene_description + " from a multicopy orthogroup " + s_orthogroup_title + ", the median coverage is " + str(d_gene_description_to_the_median_coverage_in_its_exons[s_gene_description]) + ". It was calculated using " + str(len(dl_gene_description_to_the_list_of_coverages_in_its_exons[s_gene_description])) + " positions.\n")
			
			#если ни одного однокопийного гена найдено не было (крайне маловероятно, но, в принципе, такое может быть), но многокопийные гены были, то, для простоты, считаю медианным покрытием однокопийных генов медианное покрытие по всем многокопийным генам, делённое пополам.
			if len(l_coverages_in_exons_of_single_copy_genes) == 0:
				n_median_coverage_of_exons_of_single_copy_genes = statistics.median(l_coverages_in_exons_of_genes) / 2
				f_logs.write("\nWarning! No single-copy orthogroups were found. Hence, as the approximate coverage of genes in single-copy orthogroups I take half the median coverage by positions of genes from multicopy, which is " + str(n_median_coverage_of_exons_of_single_copy_genes) + "\n\n")
			#если хотя бы один однокопийный ген был.
			else:
				n_median_coverage_of_exons_of_single_copy_genes = statistics.median(l_coverages_in_exons_of_single_copy_genes)
				f_logs.write("\nThe median coverage in exons of genes from single-copy BUSCO orthogroups is " + str(n_median_coverage_of_exons_of_single_copy_genes) + ". It was calculated using " + str(len(l_coverages_in_exons_of_single_copy_genes)) + " positions.\n\n")

			#теперь иду по всем ортогруппам и для каждой многокопийной ортогруппы считаю среднее покрытие в ней. Если оно < 0.75*(медианное покрытие в однокопийных генах), то считаю, что это ложно многокопийных ортогруппа. А если >=0.75*(медианное покрытие в однокопийных генах), то считаю, что истинно многокопийная. Параллельно, считаю количество генов в истинно многокопийных ортогруппах и ложно многокопийных ортогруппах.
			n_number_of_true_multicopy_genes = 0
			n_number_of_false_multicopy_genes = 0
			
			for s_orthogroup_title in dl_orthogroup_title_to_the_list_of_its_genes:
				#если это многокопийная ортогруппа
				if len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title]) > 1:
					n_mean_coverage_of_genes_in_this_orthogroup = 0
					for s_gene_description in dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title]:
						n_mean_coverage_of_genes_in_this_orthogroup += d_gene_description_to_the_median_coverage_in_its_exons[s_gene_description] / len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])
					
					if n_mean_coverage_of_genes_in_this_orthogroup < 0.75 * n_median_coverage_of_exons_of_single_copy_genes:
						f_logs.write("The mean coverage of genes from a multicopy BUSCO orthogroup " + s_orthogroup_title + " which contains " + str(len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])) + " genes is " + str(round(n_mean_coverage_of_genes_in_this_orthogroup, 1)) + ". It is smaller than " + str(round(0.75 * n_median_coverage_of_exons_of_single_copy_genes, 1)) + ", hence this orthogroup is considered a false multicopy.\n")
					
						n_number_of_false_multicopy_genes += len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])
						
					#в принципе, это условие и следующее можно объединить в одно (">="). Но для удобства чтения логов я разделю случай ">" и случай "=". Впрочем, думаю, случай "=" будет встречаться крайне редко.
					elif n_mean_coverage_of_genes_in_this_orthogroup == 0.75 * n_median_coverage_of_exons_of_single_copy_genes:
						f_logs.write("The mean coverage of genes from a multicopy BUSCO orthogroup " + s_orthogroup_title + " which contains " + str(len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])) + " genes is " + str(round(n_mean_coverage_of_genes_in_this_orthogroup, 1)) + ". It is equal to 0.75 * (median_coverage_in_exons_of_single_copy_genes), hence this orthogroup is considered a true multicopy.\n")
						
						n_number_of_true_multicopy_genes += len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])
					
					elif n_mean_coverage_of_genes_in_this_orthogroup > 0.75 * n_median_coverage_of_exons_of_single_copy_genes:
						f_logs.write("The mean coverage of genes from a multicopy BUSCO orthogroup " + s_orthogroup_title + " which contains " + str(len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])) + " genes is " + str(round(n_mean_coverage_of_genes_in_this_orthogroup, 1)) + ". It is larger than " + str(round(0.75 * n_median_coverage_of_exons_of_single_copy_genes, 1)) + ", hence this orthogroup is considered a true multicopy.\n")
						
						n_number_of_true_multicopy_genes += len(dl_orthogroup_title_to_the_list_of_its_genes[s_orthogroup_title])
			
			n_AG = n_number_of_single_copy_genes_found_in_the_assembly + n_number_of_true_multicopy_genes
	else:
		f_logs.write("AG is 0. Number of genes in single-copy orthogroups is 0. Number of genes in true multicopy orthogroups is 0. Number of genes in false multicopy orthogroups is 0.\n")
		f_AG_calculation_results.write("AG is 0")
		sys.exit()

	f_logs.write("AG is " + str(n_AG) + ". Number of genes in single-copy orthogroups is " + str(n_number_of_single_copy_genes_found_in_the_assembly) + ". Number of genes in true multicopy orthogroups is " + str(n_number_of_true_multicopy_genes) + ". Number of genes in false multicopy orthogroups is " + str(n_number_of_false_multicopy_genes) +".\n")
	f_AG_calculation_results.write("AG is " + str(n_AG))
	
	f_logs.close
	
	#Строю синаплот с покрытием генов.
	os.system("python3 " + s_path_to_the_folder_where_Mabs_lies + "/Additional/plot_gene_coverage_distribution.py " + s_path_to_the_output_folder + "/logs.txt 2.5 auto " + s_path_to_the_output_folder + "/gene_coverage_distribution")


