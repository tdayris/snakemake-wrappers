from irods.column import Criterion
from irods.models import DataObject, DataObjectMeta, Collection, CollectionMeta
from irods.session import iRODSSession
import ssl
import os
import pandas
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Let user provide specific irods env file
env_path = snakemake.input.get("env_file", '~/.irods/irods_environment.json')
env_file = os.path.expanduser(env_path)

logging.info("Env file found at: %s", env_file)

ssl_context = ssl.create_default_context(
    purpose=ssl.Purpose.SERVER_AUTH,
    cafile=None,
    capath=None,
    cadata=None
)
ssl_settings = {'ssl_context': ssl_context}


# Let user input multiple source of samples
sample_list = snakemake.params.get("samples", [])
if (samples_file := snakemake.input.get("samples", None)) is not None:
    with open(samples_file, "r") as samples_stream:
        for sample in samples_stream:
            sample_list.append(sample[:-1])
sample_list = list(set(sample_list))  # Drop duplicates
logging.debug("Sample list defined as: %s", str(sample_list))

# Open session
with iRODSSession(irods_env_file=env_file, **ssl_settings) as session:
    query = session.query(
        Collection.name,      DataObject.id,       DataObject.name,
        DataObject.size,      DataObject.checksum, DataObjectMeta
    )
    logging.info("Query performed")
    result_dict = {}

    # Format query result
    for res in query:
        key = (res[DataObject.id])
        try:
            result_dict[key][str(res[DataObjectMeta.name])] = str(res[DataObjectMeta.value])
        except KeyError:
            result_dict[key] = {
                "path": "{}/{}".format(res[Collection.name], res[DataObject.name]),
                "size": res[DataObject.size],
                "checksum": res[DataObject.checksum],
                str(res[DataObjectMeta.name]): str(res[DataObjectMeta.value])
            }

    # Filter query
    results = pandas.DataFrame.from_dict(result_dict, orient='index')
    logging.info("Query formatted")
    logging.debug(results.head())
    results.reset_index(inplace=True)

    results.to_csv(snakemake.output[0], sep="\t")