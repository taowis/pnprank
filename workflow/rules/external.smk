# Inputs that come from outside: reference files and the atlas downloads.
# The truth-set source snapshots (HGNC, PanelApp, ClinGen) are made by
# `pnprank build-candidates` during Stage 1 and are not rebuilt here, because that
# command also rewrites the curated candidate list.
import yaml

RAW = config["data_root"] + "/raw"
EXT = config["external"]
ATLAS_REGISTRY = yaml.safe_load(open("config/atlases.yaml"))["atlases"]


def registered_url(atlas, name):
    a = next(x for x in ATLAS_REGISTRY if x["id"] == atlas)
    f = next(x for x in a.get("files") or [] if x["name"] == name)
    return f["url"]


def atlas_raw_files(atlas):
    """Raw files of one atlas that the pipeline was run on (registry: downloaded + url)."""
    a = next(x for x in ATLAS_REGISTRY if x["id"] == atlas)
    return [f"{RAW}/{atlas}/{f['name']}" for f in a.get("files") or []
            if f.get("downloaded") and f.get("url")]


rule fetch_external:
    output:
        EXT["gnomad_constraint"]["dest"],
        EXT["gtex_median_tpm"]["dest"],
        EXT["mouse_gtf"]["dest"],
        EXT["mouse_gtf"]["genes_table"],
        EXT["human_gtf"]["dest"],
        EXT["ucsc_cytoband"]["dest"],
        config["reference"]["orthology"]["dest"],
        config["reference"]["orthology"]["table"],
    shell:
        "python workflow/scripts/fetch_external.py"


rule download_atlas_file:
    output:
        RAW + "/{atlas}/{name}",
    wildcard_constraints:
        atlas="[^/]+",
    params:
        url=lambda w: registered_url(w.atlas, w.name),
    shell:
        "python workflow/scripts/download_atlas.py {wildcards.atlas} {wildcards.name}"
