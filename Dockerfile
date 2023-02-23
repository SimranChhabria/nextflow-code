FROM ubuntu:22.04
COPY . .
ENV NEXTFLOWOUT /workspace
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
ENV PATH ./nextflow:$PATH
#CMD ["nextflow", "/workspace/nf-oncology/main.nf", "--run_name", "$RUN_NAME", "--outDir", "$NEXTFLOWOUT"]
CMD ["which","nextflow"]