use std::cmp::Ordering;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{bail, Context};
use burn::data::dataloader::batcher::Batcher;
use burn::module::Module;
use burn::prelude::Backend;
use burn::record::{CompactRecorder, Recorder};
use burn::tensor::activation::softmax;
use clap::{Args, Subcommand};
use itertools::Itertools;
#[cfg(any(feature = "tch", feature = "candle"))]
use log::warn;
use log::{debug, error, info};
use rust_htslib::bam::{self, Read};

use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
use mod_kit::util::{get_ticker, Region, TAB};

use crate::batcher::PileupBatcher;
use crate::feature_loader::FeatureLoader;
use crate::features::CountsFeatures;
use crate::models::conv_lstm::{ConvLstmModel, ConvLstmModelConfig};
use crate::prob_merger::ProbMerger;
use crate::range_feeder::RangeFeeder;
use crate::util::{filenames, InferenceMessage, ModelConfiguration};

#[derive(Subcommand)]
pub enum OpenChromatin {
    #[command(name = "predict")]
    Predict(EntryInference),
}

impl OpenChromatin {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            OpenChromatin::Predict(x) => x.run(),
        }
    }
}

#[derive(Args, Copy, Clone, Debug)]
pub struct BatchSettings {
    /// Collect this many windows of data for each run through the model.
    #[arg(long, default_value_t = 1024, hide_short_help = true)]
    pub batch_size: usize,
    /// Number of "batches" to collect at once, see documentation for exact
    /// calculation or run with --dryrun to see output
    #[arg(long, default_value_t = 100, hide_short_help = true)]
    pub super_batch_size: usize,
    /// Number of base pairs to step along the genome or genomic intervals to
    /// make predictions. Smaller numbers will result in better resolution but
    /// take more computation.
    #[arg(long, default_value_t = 25)]
    pub step_size: u32,
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryInference {
    /// Input modBAM with 6mA base modification calls.
    in_bam: PathBuf,
    /// Path to directory with open-chromatin model.
    #[arg(long = "model")]
    model_path: PathBuf,
    #[command(flatten)]
    batch_settings: BatchSettings,
    #[arg(long, short = 't', default_value_t = 4)]
    /// Require this many reads covering each prediction window.
    #[arg(long, default_value_t = 5)]
    min_coverage: u32,
    /// Only emit records/windows where the open chromatin probability is
    /// greater than or equal to this value.
    #[arg(long = "threshold")]
    filter_threshold: Option<f32>,
    /// Output bedGraph file, or "-"/"stdout" to write to stdout
    #[arg(long, short = 'o')]
    output: String,
    /// Force overwrite output file.
    #[arg(long, default_value_t = false)]
    force: bool,
    /// Path to file to write logs to, setting this parameter is recommended.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// GPU device index (number) or "cpu" for CPU. On Apple, use "mps" for
    /// Apple GPU.
    #[cfg(any(feature = "tch", feature = "candle"))]
    #[arg(short = 'x', long)]
    device: Option<String>,
    /// BED file of regions over which to make predictions
    #[arg(long, conflicts_with = "region")]
    include_bed: Option<PathBuf>,
    /// Region in the format <chrom>:start-end over which to make predictions
    #[arg(long, conflicts_with = "include_bed")]
    region: Option<String>,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Don't print header in output
    #[arg(long, default_value_t = false)]
    no_header: bool,
    /// Perform setup and validations, but stop before performing inference.
    #[arg(long, default_value_t = false)]
    dryrun: bool,
}

impl EntryInference {
    fn validate_model_directory(&self) -> anyhow::Result<()> {
        let model_config_fp = self.model_path.join(filenames::CONFIG_FN);
        if !model_config_fp.exists() {
            bail!("failed to find model configuration at {model_config_fp:?}")
        }
        let model_weights_fp = self.model_path.join(filenames::WEIGHTS_FN);
        if !model_weights_fp.exists() {
            bail!("failed to find model weights at {model_config_fp:?}")
        }
        Ok(())
    }
    fn load_model_config(&self) -> anyhow::Result<ModelConfiguration> {
        ModelConfiguration::from_path(
            &self.model_path.join(filenames::CONFIG_FN),
        )
    }

    fn load_model<B: Backend>(
        &self,
        model_config: &ModelConfiguration,
        device: &B::Device,
    ) -> anyhow::Result<ConvLstmModel<B>> {
        let weights_fp = self.model_path.join(filenames::WEIGHTS_FN);
        info!("loading weights from {:?}", &weights_fp);
        let model_weights = CompactRecorder::new()
            .load(weights_fp, device)
            .context("failed to load model weights")?;
        let model: ConvLstmModel<B> = ConvLstmModelConfig::new(
            model_config.num_features,
            model_config.hidden_size,
            model_config.num_classes,
        )
        .init(None, device)
        .load_record(model_weights);
        Ok(model)
    }

    pub fn run(&self) -> anyhow::Result<()> {
        #[cfg(all(feature = "tch", feature = "candle"))]
        compile_error!(
            "cannot have 'tch' and 'candle' features at the same time"
        );

        let _ = modkit_logging::init_logging(self.log_filepath.as_ref());
        let filter_threshold = match self.filter_threshold.as_ref() {
            Some(f) if *f < 0f32 => {
                bail!("--filter-threshold cannot be less than zero")
            }
            Some(f) if *f >= 1.0f32 => {
                bail!("--filter-threshold must be less than 1.0")
            }
            Some(f) => *f,
            None => 0f32,
        };
        self.validate_model_directory()?;

        #[cfg(all(feature = "tch", not(feature = "candle")))]
        type BurnBackend = burn::backend::LibTorch;
        #[cfg(all(feature = "candle", not(feature = "tch")))]
        type BurnBackend = burn::backend::Candle;
        #[cfg(not(any(feature = "tch", feature = "candle")))]
        type BurnBackend = burn::backend::ndarray::NdArray;

        #[cfg(all(feature = "tch", not(feature = "candle")))]
        let device = if let Some(raw_device) = self.device.as_ref() {
            match raw_device.as_str() {
                "mps" | "MPS" => burn::backend::libtorch::LibTorchDevice::Mps,
                "cpu" | "CPU" => burn::backend::libtorch::LibTorchDevice::Cpu,
                _ => {
                    if let Ok(index) = raw_device.parse::<usize>() {
                        burn::backend::libtorch::LibTorchDevice::Cuda(index)
                    } else {
                        bail!(
                            "invalid device {raw_device}, should be a number \
                             for CUDA or 'MPS' for apple"
                        )
                    }
                }
            }
        } else {
            warn!("no device given, but gpu features enabled");
            burn::backend::libtorch::LibTorchDevice::default()
        };

        #[cfg(all(feature = "candle", not(feature = "tch")))]
        let device = if let Some(raw_device) = self.device.as_ref() {
            match raw_device.as_str() {
                "mps" | "MPS" => bail!("MPS not available for Candle backend"),
                "cpu" | "CPU" => {
                    warn!("Candle CPU performance may not be ideal");
                    burn::backend::candle::CandleDevice::default()
                }
                _ => {
                    if let Ok(index) = raw_device.parse::<usize>() {
                        burn::backend::candle::CandleDevice::cuda(index)
                    } else {
                        bail!(
                            "invalid device {raw_device}, should be a number."
                        )
                    }
                }
            }
        } else {
            warn!("Candle CPU performance may not be ideal");
            burn::backend::candle::CandleDevice::default()
        };

        #[cfg(not(any(feature = "tch", feature = "candle")))]
        let device = burn::backend::ndarray::NdArrayDevice::Cpu;

        info!("using device {device:?}");

        let mpb = indicatif::MultiProgress::new();
        let reader = bam::Reader::from_path(&self.in_bam)?;
        let header = reader.header();
        let region = self
            .region
            .as_ref()
            .map(|r| Region::parse_str(r, &header))
            .transpose()?;
        if let Some(r) = region.as_ref() {
            info!("running inference on {r}");
        }

        let mut writer: Box<dyn Write + Send> = {
            if (self.output.as_str() == "-")
                || (self.output.as_str() == "stdout")
            {
                Box::new(BufWriter::new(std::io::stdout()))
            } else {
                let fp = std::path::Path::new(&self.output);
                if self.force {
                    if fp.exists() {
                        info!("overwriting output file {fp:?}");
                    }
                    Box::new(BufWriter::new(
                        std::fs::File::create(fp)
                            .context(format!("failed to create {fp:?}"))?,
                    ))
                } else {
                    Box::new(BufWriter::new(
                        std::fs::File::create_new(fp)
                            .context(format!("failed to create {fp:?}"))?,
                    ))
                }
            }
        };
        if !self.no_header {
            let header = format!("#chrom{TAB}start{TAB}stop{TAB}prob\n");
            writer.write(header.as_bytes())?;
        }
        let model_config = self.load_model_config()?;
        let model_config_str = serde_json::to_string_pretty(&model_config)?;
        mpb.suspend(|| {
            info!("loaded model config\n{model_config_str}");
        });
        let model = self.load_model::<BurnBackend>(&model_config, &device)?;
        // todo remove.
        let caller = MultipleThresholdModCaller::new_passthrough();

        let range_feeder = RangeFeeder::new(
            &self.in_bam,
            region.as_ref(),
            self.include_bed.as_ref(),
            model_config.chunk_size,
            self.batch_settings.step_size,
            self.batch_settings.batch_size as u32,
            self.batch_settings.super_batch_size as u32,
        )?;
        let feature_loader = FeatureLoader::new(
            &self.in_bam,
            caller,
            &model_config,
            model_config.chunk_size,
            self.batch_settings.step_size,
            self.min_coverage,
        )?;
        let batcher = PileupBatcher::<BurnBackend>::new(device.clone());

        if self.dryrun {
            mpb.suspend(|| info!("dryrun true - exiting."));
            std::process::exit(0);
        }

        // step 1: get batches of genome ranges that can all asynchronously be
        // run through the model
        let (snd_reg, rcv_reg) = crossbeam::channel::bounded(10);
        let regions_handle = std::thread::spawn(move || {
            for message in range_feeder {
                snd_reg.send(message).unwrap()
            }
            snd_reg.send(InferenceMessage::Done).unwrap();
        });

        // step 2: collect features for each region in the batch
        let (snd_counts, rcv_counts) = crossbeam::channel::bounded(10);
        let features_handle = std::thread::spawn(move || {
            feature_loader.run(rcv_reg, snd_counts);
        });

        // step 3: load the features onto the GPU
        let (snd_features, rcv_features) = crossbeam::channel::bounded(10);
        let batcher_device = device.clone();
        let batcher_handle = std::thread::spawn(move || {
            loop {
                match rcv_counts.recv() {
                    Ok(InferenceMessage::Work(counts)) => {
                        // debug!("{counts:?}")
                        if !counts.is_empty() {
                            let model_input =
                                batcher.batch(counts.clone(), &batcher_device);
                            snd_features
                                .send(InferenceMessage::Work((
                                    model_input,
                                    counts,
                                )))
                                .unwrap();
                        }
                    }
                    Ok(InferenceMessage::Done) => {
                        // debug!("batcher received Done, sending it on..");
                        snd_features.send(InferenceMessage::Done).unwrap()
                    }
                    Err(e) => {
                        debug!("batcher received {e}, ending");
                        break;
                    }
                }
            }
            snd_features.send(InferenceMessage::Done).unwrap();
            debug!("finishing batch loop");
        });

        // step 5 (needs to be spawned before step 4), collect the probabilities
        // and the features and write them down
        let (snd_infer, rcv_infer) = crossbeam::channel::bounded::<
            InferenceMessage<(Vec<f32>, Vec<CountsFeatures<()>>)>,
        >(10);
        let written_pb = mpb.add(get_ticker());
        written_pb.set_message("records written");
        let write_handle = std::thread::spawn(move || {
            let mut acc = ProbMerger::new(filter_threshold);
            loop {
                match rcv_infer.recv() {
                    Ok(InferenceMessage::Work((y_probs, counts))) => {
                        acc.add_preds(&counts, &y_probs);
                    }
                    Ok(InferenceMessage::Done) => {
                        debug!(
                            "writer got Done message, packing up {} records",
                            acc.size()
                        );
                        let finished = std::mem::replace(
                            &mut acc,
                            ProbMerger::new(filter_threshold),
                        );
                        let finished = finished.into_records();

                        let n = finished.len() as u64;
                        finished
                            .into_iter()
                            .sorted_by(|a, b| match a.chrom.cmp(&b.chrom) {
                                Ordering::Equal => a.start.cmp(&b.start),
                                o @ _ => o,
                            })
                            .for_each(|rec| {
                                match writer
                                    .write(format!("{rec}\n").as_bytes())
                                {
                                    Ok(_) => {}
                                    Err(e) => {
                                        error!("failed to write!, {e}")
                                    }
                                }
                            });
                        written_pb.inc(n);
                    }
                    Err(_e) => {
                        // info!("write handle got {e}");
                        break;
                    }
                }
            }
        });

        // step 4: inference, get the features loaded on the GPU, and get
        // probabilities of open chrom
        loop {
            match rcv_features.recv() {
                Ok(InferenceMessage::Work((model_input, counts))) => {
                    let [batch_size, _, _] = model_input.samples.dims();
                    let samples = model_input.samples.unsqueeze_dim(1);
                    let y_logit = model.forward(samples);
                    let y_proba = softmax(y_logit.clone(), 1);
                    let pos_proba =
                        y_proba.slice([0..batch_size, 1..2]).squeeze::<1>(1);
                    let y_probs = pos_proba.to_data().to_vec().unwrap();
                    snd_infer
                        .send(InferenceMessage::Work((y_probs, counts)))
                        .unwrap()
                }
                Ok(InferenceMessage::Done) => {
                    // debug!("inference got done message, passing it by..");
                    snd_infer.send(InferenceMessage::Done).unwrap()
                }
                Err(_e) => {
                    info!("model received completion signal");
                    break;
                }
            }
        }
        drop(snd_infer);

        regions_handle.join().unwrap();
        features_handle.join().unwrap();
        batcher_handle.join().unwrap();
        write_handle.join().unwrap();

        Ok(())
    }
}
