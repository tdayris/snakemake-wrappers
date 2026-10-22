FROM python:3.13.16-slim
LABEL Version="0.15"

COPY .condarc /etc/conda/condarc
COPY environment_google.yaml /opt/environment_google.yaml
COPY environment_slurm.yaml /opt/environment_slurm.yaml
COPY environment_local.yaml /opt/environment_local.yaml
COPY environment_google.linux-64.pin.txt /opt/environment_google.linux-64.pin.txt
COPY environment_slurm.linux-64.pin.txt /opt/environment_local.linux-64.pin.txt
COPY environment_local.linux-64.pin.txt /opt/environment_slurm.linux-64.pin.txt
COPY Miniforge3-Linux-x86_64.sh /opt/Miniforge3-Linux-x86_64.sh
COPY apptainer_1.5.4-trixie+_amd64.deb /opt/apptainer_1.5.4-trixie+_amd64.deb
COPY apptainer-suid_1.5.4-trixie+_amd64.deb /opt/apptainer-suid_1.5.4-trixie+_amd64.deb
COPY install_dependencies_in_container.sh /opt/install_dependencies_in_container.sh
COPY init_bashrc.sh /opt/init_bashrc.sh

ENV CONDARC="/etc/conda/condarc" \
    LANG=C.UTF-8 \
    SHELL="/bin/bash"

ARG MODE="local"

RUN bash /opt/install_dependencies_in_container.sh --mode "${MODE}" --force

ENTRYPOINT ["/opt/init_bashrc.sh"]
CMD ["bash"]
