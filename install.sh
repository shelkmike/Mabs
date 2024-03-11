#! /usr/bin/env bash

<<COMMENT
This script compiles additional programs provided with Mabs.
Just run "bash install.sh" .

The folder "Additional_src" contains the source code of programs used by Mabs. When "bash install.sh" is run, the content of this folder is copied into the folder "Additional" and then compiled there. Keeping the non-compiled code in a separate folder allows for a simple re-installion of Mabs in case some problem happens during installation. To re-install Mabs, one may just run "bash install.sh" again. This command deletes the folder "Additional", then re-creates it and re-compiles everything.
COMMENT

#With this command, the script stops after the first error it encounters.
set -e

#Checking that "make" is in $PATH .
if ! command -v make &> /dev/null
then
    echo 'The path to "make" is not in $PATH . Please, add it.'
    exit
fi

#Checking that "gcc" is in $PATH .
if ! command -v gcc &> /dev/null
then
    echo 'The path to "gcc" is not in $PATH . Please, add it.'
    exit
fi

#Checking that "g++" is in $PATH .
if ! command -v g++ &> /dev/null
then
    echo 'The path to "g++" is not in $PATH . Please, add it.'
    exit
fi

#Checking that "python3" is in $PATH .
if ! command -v python3 &> /dev/null
then
    echo 'The path to "python3" is not in $PATH . Please, add it.'
    exit
fi

#Checking that "perl" is in $PATH .
if ! command -v python3 &> /dev/null
then
    echo 'The path to "perl" is not in $PATH . Please, add it.'
    exit
fi

#If the user runs this script from another folder, I'm switching to the folder where this script lies.
path_to_the_folder_with_Mabs=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $path_to_the_folder_with_Mabs

#Checking if the folder "Additional" exists already. If it does, deleting it. This allows Mabs to be re-installed by running install.sh
if [ -d "./Additional" ] 
then
    rm -rf ./Additional
fi

#Copying the content of "Additional_src" to "Additional"
cp -rp ./Additional_src ./Additional

#Changing permissions of all files to 755. If a user copied the files of Mabs to a Linux machine through a Windows machine, all files will lost their permissions, which will lead to errors.
chmod --recursive 755 ./Additional/

#Installing HMMER
cd ./Additional/HMMER
./configure
make
cd ../..

#Installing Modified_hifiasm. It is a special version of Hifiasm made for Mabs. Compared to the ordinary Hifiasm, it has an additional option "--only-primary" that forces an assembly to terminate after the GFA file with primary contigs has been made. This preliminary termination allows to save some time if only primary contigs are needed.
cd ./Additional/Modified_hifiasm
make
cd ../..

#Installing SeqTk
cd ./Additional/SeqTk
make
cd ../..

#Installing Flye
cd ./Additional/Flye
make
cd ../..

#Making mabs-hifiasm.py, mabs-flye.py and calculate_AG.py executable
chmod 755 mabs-hifiasm.py mabs-flye.py calculate_AG.py