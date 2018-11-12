FROM dynverse/dynwrap:py3.6

RUN pip install python-igraph louvain # igraph and louvain do not get installed by scanpy

RUN pip install scanpy

RUN pip install fa2

LABEL version 0.1.8

ADD . /code

ENTRYPOINT python /code/run.py
