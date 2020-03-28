From conda/miniconda3-centos7

WORKDIR /app
COPY . /app

RUN yum -y update
RUN yum install -y libXext libSM libXrender

RUN conda update -n base -c defaults conda
RUN conda install -c rdkit rdkit -y
RUN conda install -c rdkit rdkit-postgresql -y
RUN pip install -r /app/requirements.txt

RUN mkdir -p /tmp/result-mols
RUN mkdir -p /mnt/pgdata

RUN adduser zilliz
RUN yum install -y wget
RUN wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/1.7/gosu-amd64" \
    && chmod +x /usr/local/bin/gosu \
    && gosu zilliz true
RUN initdb -D /mnt/pgdata
RUN cd /mnt/pgdata
RUN pg_ctl -D /mnt/pgdata -l rdclog start

USER root

RUN cd /app
CMD python src/app.py