FROM python:3.13.16-slim AS spython-base
ADD .condarc    /etc/conda/condarc
ADD environment.yaml    /opt/environment.yaml
ADD environment.linux-64.pin.txt    /opt/environment.linux-64.pin.txt
ADD Miniforge3-Linux-x86_64.sh  /opt/Miniforge3-Linux-x86_64.sh
ADD apptainer_1.5.4-trixie+_amd64.deb   /opt/apptainer_1.5.4-trixie+_amd64.deb
ADD apptainer-suid_1.5.4-trixie+_amd64.deb  /opt/apptainer-suid_1.5.4-trixie+_amd64.deb
LABEL Version v0.12
ENV CONDARC="/etc/conda/condarc"
ENV PATH="/root/miniforge3/bin:${PWD}/.snakemake/pip-deployments/bin:${PATH}"
ENV PYTHONPATH="${PYTHONPATH:-}:${PWD}/.snakemake/pip-deployments"
ENV LANG=C.UTF-8
ENV SHELL="/bin/bash"
RUN set -euo 'pipefail'
RUN apt update
RUN apt upgrade -y
RUN apt install -y build-essential wget bash bzip2 ca-certificates gnupg2 squashfs-tools git
RUN apt install -y "/opt/apptainer_1.5.4-trixie+_amd64.deb" "/opt/apptainer-suid_1.5.4-trixie+_amd64.deb"
RUN apt install --fix-broken
RUN apt clean
RUN apt autoremove -y
RUN bash "/opt/Miniforge3-Linux-x86_64.sh" -bfc
RUN eval "$(/root/miniforge3/bin/conda 'shell.bash' 'hook')"
RUN conda create -n "snakemake" --file /opt/environment.linux-64.pin.txt
RUN conda clean --all --yes
RUN export CONDARC="/etc/conda/condarc"
RUN python3 -m pip install --root-user-action 'ignore' \
--target "${PWD}/.snakemake/pip-deployments" \
'snakemake' \
'pandas==2.3.3' \
'snakemake-storage-plugin-gcs==1.1.4' \
'snakemake-storage-plugin-http==0.3.1' \
'snakemake-storage-plugin-fs==1.1.3' \
'snakemake-software-deployment-plugin-conda==0.8.4' \
'snakemake-software-deployment-plugin-container==0.7.5' \
'snakemake-executor-plugin-slurm==2.8.0' \
'snakemake-executor-plugin-googlebatch==0.5.1'
RUN echo 'eval "$(/root/miniforge3/bin/conda shell.bash hook)"' >> "/etc/bash.bashrc"
RUN echo '&& conda activate snakemake' >> "/etc/bas.bashrc"
RUN sed -i 's|^\[|#\[|' "/etc/bash.bashrc"
RUN cat "/etc/bash.bashrc"
RUN bash "/etc/bash.bashrc"
RUN python3 -m pip cache purge
RUN ln -sfv "/root/miniforge3/bin/conda" "/usr/local/bin/conda"
RUN ln -sfv "/root/miniforge3/bin/mamba" "/usr/local/bin/mamba"
RUN rm -fv "/opt/apptainer_1.5.4-trixie+_amd64.deb" \
"/opt/apptainer-suid_1.5.4-trixie+_amd64.deb" \
"/opt/Miniforge3-Linux-x86_64.sh"
RUN echo "set -euo 'pipefail'" >> /entrypoint.sh
RUN echo "echo "snakemake-conda-apptainer container: "" >> /entrypoint.sh
RUN echo "date" >> /entrypoint.sh
RUN echo "bash /etc/bash.bashrc" >> /entrypoint.sh
RUN echo "exec "$@"" >> /entrypoint.sh
RUN echo "# Ensure all commands will stop on any error" >> /entrypoint.sh
RUN echo "# Ensure conda/snakemake are available in PATH" >> /entrypoint.sh
RUN echo "# Test conda/apptainer/snakemake availability" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
RUN echo "set -euo 'pipefail'" >> /tests.sh
RUN echo "bash /etc/bash.bashrc" >> /tests.sh
RUN echo "echo "${CONDARC}"" >> /tests.sh
RUN echo "conda config --show-sources" >> /tests.sh
RUN echo "snakemake --help" >> /tests.sh
RUN echo "snakemake --version" >> /tests.sh
RUN echo "conda --version" >> /tests.sh
RUN echo "apptainer --version" >> /tests.sh
RUN echo "pip list" >> /tests.sh
RUN echo "which conda" >> /tests.sh
RUN echo "which apptainer" >> /tests.sh
RUN echo "which snakemake" >> /tests.sh
RUN echo "which bash" >> /tests.sh
RUN echo "ls -lrth /usr/bin/bash" >> /tests.sh
RUN echo "ls -lrth /bin/bash" >> /tests.sh
RUN echo "which sh" >> /tests.sh
RUN echo "ls -lrth /bin/sh" >> /tests.sh
RUN chmod u+x /tests.sh
CMD /bin/bash /entrypoint.sh
