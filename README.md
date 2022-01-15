# terogen
Implementation of the terpenoids generation used in "Bio-inspired Chemical Space Exploration of Terpenoids"
=======================================================================================================================================

#### Introduction
------------
Source code used to create, train and sample a bio-inspired terpenoids generation model described in [Bio-inspired Chemical Space Exploration of Terpenoids](https://doi.org/10.26434/chemrxiv-2022-0l482).

The model is divided into two relatively independent parts, a Reactor, which conducts metadynamics simulations to explore the reaction space of given carbocations and a Decorator, which predicts the decorating sites and groups for a given skeleton.


#### Set up
------------
For Reactor [xtb](https://xtb-docs.readthedocs.io/en/latest/contents.html) is required. In some scripts of Decorator Spark 2.4 is required (and thus Java 8).
Reactor was tested on Linux with CPUs and Decorator was tested with a Tesla V-100.
You will wish to install TeroGen in a virtual environment to prevent conflicting dependencies.

```python

conda create -n terogen python==3.6
conda activate terogen
sh install.sh
```
#### General Usage
-------------

### Carbocation Reactor (`./reactor`)

Any arbitrary molecule can be used as the initial structure for metadynamics sumilations, while herein for tergenoids generation, the isoprenoid carbocations were used.

`reactor.sh`: This script will conduct medadynamics sumilations with specific initial structure and parameters. The output is a reactant-product list with energetics properties in tsv format. It will take about **** on 12 CPUs (Xeon E5-2609 1.70GHz) to run the demo provided in the script.
`deprotonation.py`: This will quench carbocations by exhaustive deprotonation and also output a reactant-product list in tsv format.

It will take about 10 hours on 12 CPUs (Xeon E5-2609 1.70GHz) to run the demo provided in the script.

### Skeleton Decoration (`./decorator`)
There are two step for skeleton decoration, sites prediction and groups prediction. First, the decorating sites were predicted with Transformer model trained using [OpenNMT](https://opennmt.net/OpenNMT-py/) and [PyTorch](https://pytorch.org/) and then, teh R-groups were predicted with the RNN-based model proposed by [Arús-Pous et al.](https://github.com/undeadpixel/reinvent-scaffold-decorator).

Sites prediction:
`skeleton_extraction.py`: This is used to extract the carbon skeleton from the terpenoids structure.
`site_prediction.sh`: This script is used to train and test the sites prediction model.

Groups prediction:
### This model was analogous to the scaffold decorator proposed in "[SMILES-based deep generative scaffold decorator for de-novo drug design](https://doi.org/10.1186/s13321-020-00441-8)"
`group_prediction.sh`: This script is used to train and test the sites prediction model.
