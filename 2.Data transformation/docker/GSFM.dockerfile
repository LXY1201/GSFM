FROM jupyter/datascience-notebook:r-4.1.1

RUN conda install -c conda-forge --quiet --yes \
    'r=4.1' 'r-biocmanager=1.30.16' \
    && \
    R -e 'BiocManager::install(version = "3.14")' && \
    R -e 'BiocManager::install(c("ConsensusClusterPlus", "pheatmap", "GSVA", "GSEABase", "BiocParallel"))' && \
    pip install plotly==4.5.0 && \
    pip install rpy2 && \
    pip install cyjupyter 

RUN chmod -R go+rw /home/jovyan/ && chmod -R go+x /home/jovyan/.jupyter

WORKDIR /project
CMD [ "/bin/bash" ]
