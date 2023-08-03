FROM mapostolides/metafusion

ENV PATH="${PATH}:/MetaFusion/scripts"
RUN git clone https://github.com/pintoa1-mskcc/MetaFusion.git && \
    cd /MetaFusion/ && \
    git pull
RUN chmod +x /MetaFusion/scripts