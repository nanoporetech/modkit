use anyhow::{anyhow, bail};
use clap::Args;
use itertools::Itertools;
use linear_map::LinearMap;
use log::debug;
use rust_htslib::bam;

use crate::mod_base_code::{MOD_CODE_TO_DNA_BASE, RNA_LONG_NAME_TO_CODE};
use crate::writers::bedmethyl_header;

#[derive(Args, Clone, Debug)]
pub(crate) struct BedRModArgs {
    /// Output BedRModV2 header counts based on V2 Specification.
    /// Details can be found at https://github.com/dieterich-lab/euf-specs/blob/main/bedRModv2.pdf
    #[clap(help_heading = "Output Options")]
    #[arg(long, conflicts_with_all = ["with_header", "partition_tag"], hide_short_help = true)]
    bedrmod: bool,
    /// NCBI Taxonomic identifier, details at: https://doi.org/10.1093/database/baaa062
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("unknown"), hide_short_help = true)]
    organism: String,
    /// A valid RNA type
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("RNA"), hide_short_help = true)]
    modification_type: String,
    /// Genome assembly
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("unknown"), hide_short_help = true)]
    assembly: String,
    /// Annotation source
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("unknown"), hide_short_help = true)]
    annotation_source: String,
    /// Annotation version
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("unknown"), hide_short_help = true)]
    annotation_version: String,
    /// Sequencing platform
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", default_value_t = String::from("ont"), hide_short_help = true)]
    sequencing_platform: String,
    /// Basecalling model, override model set in BAM header.
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", hide_short_help = true)]
    basecalling: Option<String>,
    /// Link to bioinformatics workflow; program name, version, and/or call;
    /// information relevant to score, coverage, or frequency calculation; etc.
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", hide_short_help = true)]
    bioinformatics_workflow: Option<String>,
    /// Information about or link to experimental protocol and design
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", hide_short_help = true)]
    experiment: Option<String>,
    /// Databank:ID of data
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, requires = "bedrmod", hide_short_help = true)]
    external_source: Option<String>,
}

impl BedRModArgs {
    pub(crate) fn enabled(&self) -> bool {
        self.bedrmod
    }

    pub(crate) fn header(
        &self,
        bam_header: &bam::HeaderView,
    ) -> anyhow::Result<String> {
        let header = bam::Header::from_template(bam_header);
        let tags = header.to_hashmap();
        let ds_record =
            tags.get("RG").and_then(|rg_tag| find_ds_tag(rg_tag)).ok_or_else(
                || anyhow!("failed to find DS in RG record mapping"),
            )?;
        let basecaller_model = ds_record
            .split_whitespace()
            .find(|x| x.starts_with("basecall_model="))
            .map(|x| x.replace("basecall_model=", ""));
        let raw_base_mod_names = parse_ds_record(ds_record)?;
        let modification_names_info = raw_base_mod_names
            .into_iter()
            .map(|name| {
                let mod_code =
                    RNA_LONG_NAME_TO_CODE.get(name.as_str()).unwrap();
                let dna_base = MOD_CODE_TO_DNA_BASE.get(mod_code).unwrap();
                format!("{name}:{mod_code}:{dna_base}")
            })
            .join(",");

        let make_header_line =
            |k: &str, v: &str| -> String { format!("#{k}={v}") };

        let empty = "".to_string();
        let header_lines = vec![
            make_header_line("fileformat", "bedRModv2"),
            make_header_line("organism", self.organism.as_str()),
            make_header_line(
                "modification_type",
                self.modification_type.as_str(),
            ),
            make_header_line(
                "modification_names",
                modification_names_info.as_str(),
            ),
            make_header_line("assembly", self.assembly.as_str()),
            make_header_line(
                "annotation_version",
                self.annotation_version.as_str(),
            ),
            make_header_line(
                "sequencing_platform",
                self.sequencing_platform.as_str(),
            ),
            make_header_line(
                "basecalling",
                self.basecalling
                    .as_ref()
                    .or(basecaller_model.as_ref())
                    .unwrap_or_else(|| &empty)
                    .as_str(),
            ),
            make_header_line(
                "bioinformatics_workflow",
                self.bioinformatics_workflow
                    .as_ref()
                    .unwrap_or_else(|| &empty)
                    .as_str(),
            ),
            make_header_line(
                "experiment",
                self.experiment.as_ref().unwrap_or_else(|| &empty).as_str(),
            ),
            make_header_line(
                "external_source",
                self.external_source
                    .as_ref()
                    .unwrap_or_else(|| &empty)
                    .as_str(),
            ),
            bedmethyl_header(),
        ];
        let header = header_lines.join("\n");
        Ok(header)
    }
}

fn parse_model_name(raw: &str) -> anyhow::Result<Vec<String>> {
    let mod_names = raw
        .split("@")
        .nth(1)
        .map(|version_codes| {
            version_codes
                .split("_")
                .skip(1)
                .map(|x| String::from(x))
                .filter(|x| RNA_LONG_NAME_TO_CODE.contains_key(x.as_str()))
                .collect::<Vec<String>>()
        })
        .ok_or_else(|| anyhow!("invalid modbase model name {raw}"))?;
    if mod_names.is_empty() {
        bail!("failed to parse mod names from {raw}")
    } else {
        debug!("modbase model {raw} has mod-names {mod_names:?}");
        Ok(mod_names)
    }
}

fn parse_ds_record(ds_record: &str) -> anyhow::Result<Vec<String>> {
    let modbase_models_list = ds_record
        .split_whitespace()
        .find(|part| part.starts_with("modbase_models="))
        .map(|x| x.replace("modbase_models=", ""))
        .ok_or_else(|| {
            anyhow!("failed to find modbase_models in DS record {ds_record}")
        })?;
    let mod_names = modbase_models_list
        .split(",")
        .map(|model_name| parse_model_name(model_name))
        .collect::<anyhow::Result<Vec<Vec<String>>>>()?;
    Ok(mod_names.into_iter().flatten().unique().collect())
}

fn find_ds_tag(
    read_group_record: &Vec<LinearMap<String, String>>,
) -> Option<&String> {
    read_group_record
        .iter()
        .flat_map(|kvs| {
            kvs.iter().find_map(|(k, v)| if k == "DS" { Some(v) } else { None })
        })
        .nth(0)
}
