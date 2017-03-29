FROM centos:latest

RUN sh -c '/bin/echo -e "y" | yum -y update' 

RUN sh -c 'bin/echo -e "y" | yum install -y epel-release'

RUN sh -c 'bin/echo -e "y" | yum -y update'

RUN sh -c 'bin/echo -e "y" | yum install -y R git'

RUN sh -c 'bin/echo -e "y" | yum group install "Development Tools"'

RUN sh -c 'bin/echo -e "y" | yum install -y emacs'

RUN yum clean all 

RUN yum -y update

RUN yum install -y libxml2-devel.x86_64 libcurl-devel 

RUN yum install -y libpng-devel.x86_64

RUN git clone https://github.com/gosuzombie/DNF

RUN Rscript ./DNF/RCode/dependencyInstaller.R







