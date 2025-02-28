import pandas

bamdict = {}
for bam in bams:
    fp = bam.split("/")[-1]
    if "baseline" in fp:
        print(f"Working on baseline bam: {fp}")
        sample = fp.split(".baseline.")[0]
        if sample in bamdict.keys():
            bamdict[sample]["baseline"] = bam
        else:
            bamdict[sample] = {"baseline": bam}
    elif "wbc" in fp:
        print(f"Working on wbc bam: {fp}")
        sample_version_manip = fp.split(".wbc.")[0]
        manip = sample_version_manip.split("_")[-1]
        version = sample_version_manip.split("_")[-2]
        sample = sample_version_manip[:-len(f"_{manip}_{version}")]
        if sample in bamdict.keys():
            if version in bamdict[sample].keys():
                if manip in bamdict[sample][version].keys():
                    print("Everything exists in bamdict.")
                    bamdict[sample][version][manip]["wbc"] = bam
                else:
                    print(f"Missing manip {manip} in version {version} of sample {sample}.")
                    bamdict[sample][version][manip] = {"wbc": bam}
            else:
                print(f"Missing version {version} in sample {sample}. Adding manip {manip} in it.")
                bamdict[sample][version] = {manip: {"wbc": bam}}
        else:
            print(f"Missing sample {sample} in bamdict. Adding version {version} and manip {manip} in it.")
            bamdict[sample] = {version: {manip: {"wbc": bam}}}
    else:
        print(f"Working on ctc bam: {fp}")
        sample_version_manip_nb = fp.split(".ctc.")[0]
        nb = sample_version_manip_nb.split("_")[-1]
        manip = sample_version_manip_nb.split("_")[-2]
        version = sample_version_manip_nb.split("_")[-3]
        sample = sample_version_manip_nb[:-len(f"_{version}_{manip}_{nb}")]
        if sample in bamdict.keys():
            if version in bamdict[sample].keys():
                if manip in bamdict[sample][version].keys():
                    if nb in bamdict[sample][version][manip].keys():
                        print(f"Duplicate nb {nb} in bamdict!")
                        raise ValueError(bamdict)
                    else:
                        print("Everything exists in bamdict.")
                        bamdict[sample][version][manip][nb] = {"ctc": bam}
                else:
                    print(f"Missing manip {manip} in version {version} of sample {sample}, number{nb}.")
                    bamdict[sample][version][manip] = {nb : {"ctc": bam}}
            else:
                print(f"Missing version {version} in sample {sample}. Adding manip {manip} in it for number {nb}.")
                bamdict[sample][version] = {manip: {nb: {"ctc": bam}}}
        else:
            print(f"Missing sample {sample} in bamdict. Adding version {version} and manip {manip} in it for number {nb}.")
            bamdict[sample] = {version: {manip: {nb: {"ctc": bam}}}}


table = []
last_wbc = None
last_baseline = None
for sample in bamdict.keys():
    print(f"Working on {sample}: {bamdict[sample]}.")
    versions = [k for k in bamdict[sample].keys() if k != "baseline"]
    for version in versions:
        print(f"Working on {version}: {bamdict[sample][version]}.")
        manips = [m for m in bamdict[sample][version].keys() if str(m).startswith("M")]
        for manip in manips:
            print(f"Working on {manip}: {bamdict[sample][version][manip]}.")
            nbs = [n for n in bamdict[sample][version][manip].keys() if n != "wbc"]
            for nb in nbs:
                wbc = bamdict[sample][version][manip].get("wbc", last_wbc)
                ctc = bamdict[sample][version][manip][nb].get("ctc", None)
                baseline = bamdict[sample].get("baseline", last_baseline)
                try:
                    table.append({"Sample_id": sample, "Version": version, "Manip": manip, "NB": nb, "CTC": ctc, "WBC": wbc, "Baseline": baseline})
                    if wbc is not None:
                        last_wbc = wbc
                    if baseline is not None:
                        last_baseline = baseline
                except KeyError:
                    print(f"Error in {sample}, {version}, {manip}, {nb}: {bamdict[sample]}.")
                    raise


design = pandas.DataFrame.from_records(data = table)
print(design.head())

