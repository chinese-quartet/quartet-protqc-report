---

![MultiReport for Quartet Proteomics](multireport-logo.png)

---
**MultiReport is a plugin based on MultiQC, providing additional tools which are
specific to Proteomics quality control of Quartet Project.**

For more information about Quartet Project, see http://chinese-quartet.org/

For more information about MultiQC, see http://multiqc.info

## Usage

To use this plugin, you need to install MultiQC and install `quartet-proteome-report`.

```shell
# Install MultiQC
pip install multiqc

# Install quartet-proteome-report
git clone https://github.com/chinese-quartet/quartet-protqc-report.git
cd quartet-proteome-report
python setup.py install
```

*The input files for the MultiReport is the output result of `qcprot` (https://github.com/chinese-quartet/qcprot)*


Then, you can get the QC report by the following actions:

```shell
# E.g., save all data into the folder `results`
multiqc ./results/

# For the results which is not belong to Quartet Proteomics pipeline, you can use the the original MultiQC
multiqc ./results/ --disable-plugin
```

## Development
If you're developing this code, you'll want to clone it locally and install
it manually instead of using `pip`:

```shell
git clone https://github.com/chinese-quartet/quartet-protqc-report.git
cd report/quartet-proteome-report
# You don't need to rerun the installation every time you make an edit (though you still do if you change anything in setup.py).
python setup.py develop
```