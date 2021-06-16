#!/bin/sh
sudo apt-get update
sudo apt-get install -y unzip
sudo apt-get install -y build-essential
sudo apt-get install -y m4
cd ..
wget https://files.gap-system.org/gap-4.11/tar.gz/gap-4.11.0.tar.gz
tar -xzf gap-4.11.0.tar.gz
cd gap-4.11.0
./configure && make
cd pkg
git clone https://github.com/neshitov/localsnf
cd localsnf
./configure && make
cd
cd ch2bt_bld
mkdir dim_6_logs
alias gap='~/gap-4.11.0/bin/gap.sh'
