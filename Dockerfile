FROM dynverse/dynwrap:py3.6

LABEL version 0.1.3

RUN pip install python-igraph louvain # igraph and louvain do not get installed by scanpy

RUN pip install feather-format

RUN pip install scanpy

RUN pip install fa2

ADD . /code

ENTRYPOINT python /code/run.py
