FROM dynverse/dynwrappy3:v0.1.0

RUN pip install python-igraph louvain # igraph and louvain do not get installed by scanpy

RUN pip install scanpy

RUN pip install fa2

ADD . /code

ENTRYPOINT python /code/run.py
