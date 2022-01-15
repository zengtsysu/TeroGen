#!/bin/sh

####
This model was analogous to the scaffold decorator proposed in "[SMILES-based deep generative scaffold decorator for de-novo drug design](https://doi.org/10.1186/s13321-020-00441-8)"

echo " Job begin time: " `date`

##### randomize smiles
python ./create_randomized_smiles.py -i data/training.smi -o data/training -n 50

##### creat model
#python create_model.py -i data/training/001.smi -o data/checkpoints/rgroup.empty -d 0.2

##### train model(GPU)
export CUDA_VISIBLE_DEVICES=0
#python train_model.py -i terpenoids/models/model.empty -o terpenoids/models/model.trained -s terpenoids/training -e 100 -b 64 --sen 10 --collect-stats-frequency 1 --collect-stats-validation-set-path terpenoids/validation/validation.smi --collect-stats-log-path terpenoids/log


##### sampling
$spark-submit --driver-memory=38g sample_scaffolds.py -m data/models/rgroup.trained.80 -i test.txt -o test.parquet -r 10 -n 10 -d single

echo " Job finish time: " `date`
