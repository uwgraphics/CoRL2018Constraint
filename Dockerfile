FROM ubuntu:bionic
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install python2.7 git 
RUN apt-get -y install python-scipy python-sympy 
RUN apt-get -y install python-pip
RUN pip install matplotlib
RUN pip install --upgrade cython
RUN mkdir /CoRL2018Constraints
COPY . /CoRL2018Constraints
WORKDIR /CoRL2018Constraints
# generating the constraint equation library
RUN python constraint_equations/geometric_constraints.py
CMD ["python", "test.py"]