Nucleosome Dynamics - docker implementation
===================
Nucleosome Dynamics docker image is a containerized implementation of 'Nucleosome Dynamics suite' particularly developed for use with the [Galaxy platform](https://dev.usegalaxy.es/).

#### Nucleosome Dynamics
Nucleosome Dynamics is a suite of R programs for **nucleosome-related analyses** based on MNase-seq experimental data. The toolkit includes the following analyses:

- nucleR: define and classify the location of nucleosomes from MNase-seq data.
- Nucleosome Dynamics: compares different MNase-seq experiments identifying variations in nucleosomes location.
- NFR: locate nucleosome-free regions
- Periodicity: predicts  nucleosome phasing at gene level
- TSS: classify transcription start sites based on the surrounding nucleosomes
- Stiffness: computes nucleosome stiffness derived from a Gaussian function fitting

#### Container software 
This docker image extends and distributes the following rellevant software:

* R libraries:
   * NucleR R library: https://github.com/nucleosome-dynamics/nucleR
   * Nucleosome Dynamics R library: https://github.com/nucleosome-dynamics/NucDyn
* R programs:
   * Nucleosome Dynamics analyses: https://github.com/nucleosome-dynamics/nucleosome_dynamics

# Build the image
The docker image for Nucleosome Dynamics can be found in the [docker hub](https://hub.docker.com/r/mmbirb/nucleosome-dynamics) as `mmbirb/nucleosome-dynamics`. However, this repository contains the data to manually rebuild it by running **docker build**.

```sh
docker build -t nucleosome-dynamics .
```
> Note that the current working directory must contain the Dockerfile file.

# Running nucleosome-dynamics container
You can either run the container manually, or via the web-based [Galaxy platform](https://dev.usegalaxy.es/). In both cases, individual nucleosome-related analyses are offered as:
- individual analyses
- workflow

> Install [docker](https://docs.docker.com/engine/installation/) for your system, if not previously done.

### Running an individual analysis
You can run manually your containers using the following command:

```sh
docker run -v /path/to/data_dir/:/path/to/data_dir mmbirb/nucleosome-dynamics  [analysis] [analysis_options]
```

In case you do not have the container stored locally, docker will download it for you. Once the execution is finished, the container will exit.

A short description of the parameters would be:
- `docker run` will run the container for you.
- `-v /path/to/data_dir/:/path/to/data_dir` will make the `data_dir` where input files are stored available to the container
- `mmbirb/nucleosome-dynamics` is the image name, which can be found in the [docker hub](https://dev.usegalaxy.es/).
- `[analysis] [analysis_options]` correspond to the analysis type you want to run and its paremeters.

Available analysis are:

| `[analysis]` | Description |
| -------- | -------- |
| readBAM	  | Read Aligned MSase-seq BAM into a RData structure (required for further process) |
| nucleR	   | Determine the positions of the nucleosomes across the genome |
| NFR		  | Short regions depleted of nucleosomes |
| txstart	  | Classify Transciption start accordint to the properties of surrounding nucleosomes	|
| periodicity  | Periodic properties of nucleosomes inside gene bodies |
| stiffness | Aparent stiffness constant foreach nucleosome obtained by fitting the coverage to a gaussian distribution |
| nucleR_stats|	 Nucleosome call statistics|
| nucDyn          | Comparison of two diferent MNase-seq experiments to nucleosome architecture local changes |
| NFR_stats|		 Nucleosome Free Regions statistics|
| txstart_stats|	 TSS and TTS statistics|
| periodicity_stats|  Statistics on Nucleosome periodicity|
| stiffness_stats|	Statistics on stiffness|
| nucDyn_stats|	   Statistics on Nucleosome Dynamics analysis|


Each analysis has its own input files and arguments. `docker run mmbirb/nucleosome-dynamics [analysis] --help` will display such information in detail. Yet, a full usage description can be found at the [Nucleosome Dynamics repository](https://github.com/nucleosome-dynamics/nucleosome_dynamics).

#### Example

Here is an example on how to load a MNase-seq reads file using `readBAM`. It takes as input a BAM file and convert it into an RData file ready to be feed to other analyses ('nucleR', 'NFR', etc.).

```sh
docker run -v $PWD/test/data/:$PWD/test/data mmbirb/nucleosome-dynamics readBAM --input $PWD/test/data/cellcycleG2_chrII.bam --output $PWD/test/data/cellcycleG2_chrII.RData --type paired
```

# Running a workflow script
You can combine different analysis tools to build your own workflow, you can do it at Galaxy platform or you can run it manually following this command:

```sh

docker run -v /path/to/data_dir/:/path/to/data_dir mmbirb/nucleosome-dynamics  run  [WF_file_path]

```
Where:
- `WF_file_path`: is a bash file containing the Rscript calls to different analyses.

#### Example

Here is an example on how to run the test workflow file `test/scripts/wf-test.sh`, a bash file sequencially calling all 'Nucleosome Dynamics' analyses.

```
docker run -v /path/to/data_dir/:/path/to/data_dir mmbirb/nucleosome-dynamics  run /path/to/data_dir/test/scripts/wf-test.sh 
```
