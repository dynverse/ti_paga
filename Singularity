#######################################################################################
## DO NOT EDIT THIS FILE! This file was automatically generated from the dockerfile. ##
## Run dynwrap:::.container_dockerfile_to_singularityrecipe() to update this file.   ##
#######################################################################################

Bootstrap: shub

From: dynverse/dynwrap:py3.6

%labels
    version 0.1.1

%post
    chmod -R a+r /code
    chmod a+x /code
    pip install python-igraph louvain # igraph and louvain do not get installed by scanpy
    pip install feather-format
    pip install scanpy
    pip install fa2

%files
    . /code

%runscript
    exec python /code/run.py
