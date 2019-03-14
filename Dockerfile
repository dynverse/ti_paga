FROM dynverse/dynwrappy:v0.1.0

# igraph and louvain do not get installed by scanpy
RUN pip install python-igraph louvain

RUN pip install scanpy

RUN pip install fa2

COPY definition.yml run.py example.sh /code/

ENTRYPOINT ["/code/run.py"]
