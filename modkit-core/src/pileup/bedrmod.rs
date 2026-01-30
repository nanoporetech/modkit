use anyhow::anyhow;
use clap::Args;
use itertools::Itertools;
use linear_map::LinearMap;
use log::warn;
use rust_htslib::bam;
use rustc_hash::FxHashMap;

use crate::mod_base_code::{
    DnaBase, ModCodeRepr, ModifiedBasesOptions, RNA_CODES_TO_MODOMICS_NAMES,
};
use crate::writers::bedrmod_bedmethyl_header;

#[derive(Args, Clone, Debug)]
pub(crate) struct BedRModArgs {
    /// Output BedRModV2 header counts based on V2 Specification.
    /// Details can be found at https://dieterich-lab.github.io/euf-specs/bedRModv2.pdf
    #[clap(help_heading = "Output Options")]
    #[arg(long, conflicts_with_all = ["with_header", "combine_strands", "combine_mods", "phased"], requires = "modified_bases", hide_short_help = true)]
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
    /// User-specified code, format: <Code> <MODOMICS Short Name> <Primary
    /// Base>. <Code> is the base modification code described in the
    /// SAMtags or a numeric ChEBI code present in the BAM file. For example 'a
    /// m6A A' would be for <Code> 'a', <MODOMICS Short Name> m6A, and
    /// <Primary Base> A.
    #[clap(help_heading = "BedRMod Options")]
    #[arg(long, action = clap::ArgAction::Append, num_args = 3, requires = "bedrmod", hide_short_help = true)]
    modomics_code: Option<Vec<String>>,
}

impl BedRModArgs {
    pub(crate) fn enabled(&self) -> bool {
        self.bedrmod
    }

    fn parse_user_modomics_codes(
        &self,
    ) -> anyhow::Result<FxHashMap<ModifiedBasesOptions, String>> {
        let mut lu = FxHashMap::default();
        for c in self.modomics_code.as_ref().unwrap_or(&Vec::new()).chunks(3) {
            assert_eq!(c.len(), 3);
            let mod_code = ModCodeRepr::parse(&c[0])?;
            let primary_base = c[2]
                .parse::<char>()
                .map_err(|e| anyhow!("illegal DNA base {e}"))
                .and_then(|c| DnaBase::parse(c).map_err(|e| e.into()))?;
            let modomics_code = c[1].to_owned();
            lu.insert(
                ModifiedBasesOptions::new(mod_code, primary_base),
                modomics_code,
            );
        }
        Ok(lu)
    }

    pub(crate) fn header(
        &self,
        bam_header: &bam::HeaderView,
        modified_base_options: &[ModifiedBasesOptions],
    ) -> anyhow::Result<String> {
        let user_modomics_codes: FxHashMap<ModifiedBasesOptions, String> =
            self.parse_user_modomics_codes()?;
        let header = bam::Header::from_template(bam_header);
        let tags = header.to_hashmap();
        let basecaller_model = match self.basecalling.as_ref() {
            Some(basecaller) => Some(basecaller.to_owned()),
            None => {
                match tags
                    .get("RG")
                    .and_then(|rg_tag| find_ds_tag(rg_tag))
                    .ok_or_else(|| {
                        anyhow!("failed to find DS in RG record mapping")
                    }) {
                    Ok(ds_record) => {
                        let modbase_models_list =
                            parse_modbase_models(&ds_record)?;
                        let basecaller_model = ds_record
                            .split_whitespace()
                            .find(|x| x.starts_with("basecall_model="))
                            .map(|x| x.replace("basecall_model=", ""))
                            .map(|x| format!("{x},{modbase_models_list}"));
                        basecaller_model
                    }
                    Err(_) => {
                        warn!(
                            "failed to parse DS record from RG header record. \
                             Basecaller not provided so will be empty"
                        );
                        None
                    }
                }
            }
        };
        let modification_names_info = modified_base_options
            .iter()
            .map(|modified_base_option| {
                user_modomics_codes
                    .get(modified_base_option)
                    .map(|x| x.to_string())
                    .or_else(|| {
                        RNA_CODES_TO_MODOMICS_NAMES
                            .get(modified_base_option)
                            .map(|x| x.to_string())
                    })
                    .ok_or_else(|| {
                        anyhow!(
                            "need to provide modomics code for \
                             {modified_base_option:?}"
                        )
                    })
            })
            .collect::<anyhow::Result<Vec<_>>>()?
            .into_iter()
            .zip(modified_base_options.iter())
            .map(|(modomics_name, mbo)| {
                let mod_code = mbo.mod_code;
                let dna_base = mbo.primary_base;
                format!("{mod_code}:{modomics_name}:{dna_base}")
            })
            .join(",");

        let make_header_line =
            |k: &str, v: &str| -> String { format!("#{k}={v}") };

        let empty = "".to_string();
        let command_line = std::env::args().collect::<Vec<String>>().join(" ");
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
                "annotation_source",
                self.annotation_source.as_str(),
            ),
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
                basecaller_model.as_ref().unwrap_or_else(|| &empty).as_str(),
            ),
            make_header_line(
                "bioinformatics_workflow",
                self.bioinformatics_workflow
                    .as_ref()
                    .unwrap_or_else(|| &command_line)
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
            bedrmod_bedmethyl_header(),
        ];
        let header = header_lines.join("\n");
        Ok(header)
    }
}

fn parse_modbase_models(ds_record: &str) -> anyhow::Result<String> {
    ds_record
        .split_whitespace()
        .find(|part| part.starts_with("modbase_models="))
        .map(|x| x.replace("modbase_models=", ""))
        .ok_or_else(|| {
            anyhow!("failed to find modbase_models in DS record {ds_record}")
        })
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
